clear all;
close all;
restoredefaultpath();
addpath(genpath('./utils'))
addpath(genpath('./libraries'))

warning('off','geopdes:nrbmultipatch'); %-> Warns for using nrbmultipatch instead of mp_geo_load with files
warning('off','MATLAB:decomposition:LoadNotSupported'); %-> Composite structure don't allow to load decomposition objects (but are usable inside of spmd-blocks)

maxNumCompThreads(15); % Save resources if working on a Cluster

%% Setup numerical tests
degrees = 1:3;
divs = [round(2.^(1:0.5:3))];
steps = [round(2.^(1:1:13))];
subs = 1;

%% Time setup
t0 = 0;
tend = 1;

%% Materials
air_regions = 1;
copper_regions = 2;

% mu_air = 1.25663753e-6;
nu_air = 1;
sigma_air = 0;
% mu_copper = 1.256629e-6;
nu_copper = 1;
sigma_copper = 1;

nu  = @(x, y, z, ind) ((ismember(ind,air_regions))*nu_air + ...
    (~ismember(ind,air_regions))*nu_copper)*ones(size(x));
sigma = @(x, y, z, ind) ((ismember(ind,air_regions))*sigma_air + ...
    (~ismember(ind,air_regions))*sigma_copper)*ones(size(x));

nu_mat = @(x, y, z, ind) cat(1,reshape(nu(x,y,z,ind),[1,size(x)]),...
                               reshape(nu(x,y,z,ind),[1,size(x)]),...
                               reshape(nu(x,y,z,ind),[1,size(x)]));
sigma_mat = @(x, y, z, ind) cat(1,reshape(sigma(x,y,z,ind),[1,size(x)]),...
                                  reshape(sigma(x,y,z,ind),[1,size(x)]),...
                                  reshape(sigma(x,y,z,ind),[1,size(x)]));

%% Analytical solution
timeDep = @(t) exp(-t(:).');
AFun = @(x, y, z, t) timeDep(t).*cat(1,reshape(cos(y).*cos(z).*sin(x),[1,size(x)]),...
                                 reshape(-2.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                                 reshape(cos(x).*cos(y).*sin(z),[1,size(x)]));
BFun = @(x, y, z, t) timeDep(t).*cat(1,reshape(-3.*cos(x).*sin(y).*sin(z),[1,size(x)]),...
                                    zeros([1,size(x)]),...
                                    reshape(3.*cos(z).*sin(x).*sin(y),[1,size(x)]));

dtAFun = @(x, y, z, t) -AFun(x, y, z, t);
EFun = @(x, y, z, t) -dtAFun(x, y, z, t);

%% Source, initial condition, boundary conditions
initFun = @(x, y, z, ind) AFun(x,y,z,0);

cc_sol = @(x, y, z, ind) cat(1,reshape(3.*cos(y).*cos(z).*sin(x),[1,size(x)]),...
                               reshape(-6.*cos(x).*cos(z).*sin(y),[1,size(x)]),...
                               reshape(3.*cos(x).*cos(y).*sin(z),[1,size(x)])).*nu_mat(x,y,z,ind);
dt_sol = @(x, y, z, ind) -initFun(x,y,z,ind).*sigma_mat(x,y,z,ind);
sourceFunS = @(x, y, z, ind) dt_sol(x, y, z, ind) + cc_sol(x, y, z, ind);

dirFunS = @(x, y, z, ind) AFun(x,y,z,0);


%% Preparation for measurements
errBL2max = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));
errACL2max = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));
errECL2max = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));

errB_acevedo = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));
errE_acevedo = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));

numIter = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));
condEst = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));
coarseSize = NaN*zeros(numel(subs),numel(degrees),numel(divs),numel(steps));

filename = 'data/astar_cube2_results.csv';
fid = fopen(filename,'w');
fprintf(fid,'subs,deg,divs,steps,pri,cond,iter,errBmax,errAmax,errEmax,errBa,errEa\n');

for k = 1:numel(subs)

    %% Setup of subdomains
    nrbarr = cnstrct_cube_nrbarr([1,1,1],subs(k)*[2,1,1]);
    nrbcell = num2cell(nrbarr)';
    clear nrbarr;
    
    % Stuff for plotting
    mid = [0;0;0];
    for i=1:numel(nrbcell)
        mid = mid + nrbeval(nrbcell{i},[0.5,0.5,0.5]);
    end
    mid = mid./numel(nrbcell);
    
    %% Boundary and interfaces (all IETI + Dirichlet)
    [interfacesIETI,boundariesIETI] = nrbmultipatch([nrbcell{:}]);
    for i=1:numel(boundariesIETI)
        dirSides(i).patch = boundariesIETI(i).patches;
        dirSides(i).side = boundariesIETI(i).faces;
    end

    %% Load geometries
    geoCell = cellfun(@(nrb) geo_load(nrb), nrbcell,'UniformOutput', false);

    for p=1:numel(degrees)
        for s=1:numel(divs)

            %% Setup for spaces and meshes
            degree = degrees(p)*[1 1 1];
            degCell = cellfun(@(nrb) degree, nrbcell, 'UniformOutput', false);
            ndiv = divs(s)*[1 1 1];
            ndivCell = cellfun(@(nrb) ndiv, nrbcell, 'UniformOutput', false);
            nquadCell = cellfun(@(deg) deg+1, degCell, 'UniformOutput', false);
            regCell = cellfun(@(deg) deg-1, degCell, 'UniformOutput', false);
    
            [mshCell,spCurlCell,spGradCell] = cellfun(@(geo,ndiv,deg,reg,nquad) setupMeshAndSpaces(geo,ndiv,deg,reg,nquad),...
                                    geoCell,ndivCell, degCell, regCell, nquadCell, 'UniformOutput', false);
            
    
            %% Build global graph for periodic and ieti conditions
            % Setup local graphs
            indeces = 1:numel(nrbcell);
            indCell = num2cell(indeces)';
            locGraphCell = cellfun(@(sp,nrb,idx) setupLocalGraph(sp,nrb,idx,dirSides,interfacesIETI,[]),...
                spCurlCell,nrbcell, indCell, 'UniformOutput',false);
    
            % Construct global numbering (for nodes and edges)
            [gnumGrad,~] = mp_interface(interfacesIETI, spGradCell);
            [gnumCurl,~] = mp_interface_hcurl(interfacesIETI, spCurlCell);
    
            % Build global graph
            gloGraph = loc2glo_graph(locGraphCell,gnumCurl,gnumGrad,0);
            T = minspantree(gloGraph,'Method','sparse','Type','forest');
            gloTree = T.Edges.IDs;
    
            %% Setup DOF subsets
            ndofCell = cellfun(@(sp) sp.ndof, spCurlCell, 'UniformOutput', false);
            ndofCumu = cumsum([0,[ndofCell{:}]]);
    
            [gDirCell,dirCell] = cellfun(@(sp,msh,idx) setupDirichletCondition(sp,msh,dirFunS,idx,dirSides),...
                spCurlCell, mshCell, indCell, 'UniformOutput', false);
    
            ietiCell = cellfun(@(sp,msh,idx) getIETIDOFs(sp,msh,idx,interfacesIETI),...
                spCurlCell, mshCell, indCell, 'UniformOutput', false);
            ietiCell = cellfun(@(ieti,dir) setdiff(ieti,dir), ietiCell, dirCell, 'UniformOutput', false);
    
            treeCell = cellfun(@(gnum) glo2locTree(gloTree, gnum)', gnumCurl, 'UniformOutput', false);
            treeCell = cellfun(@(tree,dir) setdiff(tree,dir), treeCell, dirCell, 'UniformOutput', false); % Remove Dirichlet
    
            % Primal DOFs
            priCell = cellfun(@(tree,ieti) intersect(tree,ieti), treeCell, ietiCell, 'UniformOutput', false);
            
            % Eliminated DOFs
            eliTreeCell = cell(numel(nrbcell),1);
            for i=1:numel(nrbcell)
                if ismember(i,air_regions)
                    eliTreeCell{i} = setdiff(treeCell{i},priCell{i});
                end
            end
            eliCell = cellfun(@(dir,eliTree) [dir;eliTree], dirCell, eliTreeCell, 'UniformOutput', false);

            % Remaining DOFs
            remCell = cellfun(@(sp,pri,eli) setdiff((1:sp.ndof)',union(eli,pri)), spCurlCell, priCell, eliCell, 'UniformOutput', false);
    
            % Splitting of remaining DOFs for preconditioning
            remIntCell = cellfun(@(rem,ieti) intersect(rem,ieti), remCell, ietiCell, 'UniformOutput', false);
            remVolCell = cellfun(@(rem,remInt) setdiff(rem,remInt), remCell, remIntCell, 'UniformOutput', false);

            subsets = cell(numel(treeCell),1);
            for i = 1:numel(treeCell)
                subsets{i}.eli = eliCell{i};
                subsets{i}.rem = remCell{i};
                subsets{i}.pri = priCell{i};
    
                subsets{i}.remVol = remVolCell{i};
                subsets{i}.remInt = remIntCell{i};
            end
       
            %% Setup IETI-constraints
            [B,Cnull] = setupCouplingMatTI(gnumCurl,ndofCumu(end));
    
            %% Compute multiplier subsets
            remCell = cellfun(@(set,ndof) ndof + set.rem, subsets, num2cell(ndofCumu(1:end-1))', 'UniformOutput', false);
            [remMults,~,~] = find(B(:,vertcat(remCell{:})));
            remMults = unique(remMults);
    
            priCell = cellfun(@(set,ndof) ndof + set.pri, subsets, num2cell(ndofCumu(1:end-1))', 'UniformOutput', false);
            [priMults,~,~] = find(B(:,vertcat(priCell{:})));
            priMults = unique(priMults);
    
            multNeli = union(priMults,remMults);
            [~,multPri,~] = intersect(multNeli,priMults);
            [~,multRem,~] = intersect(multNeli,remMults);
    
            for i=1:numel(subsets)
                subsets{i}.multPri = multPri;
                subsets{i}.multRem = multRem;
            end
    
            %% Construct RHS for fully eliminated constraints
            geSCell = cellfun(@(gDir,eliTree) [gDir;zeros(numel(eliTree),1)], gDirCell, eliTreeCell, 'UniformOutput', false);
    
            %% Assemble nullspace representation
            B = B(multNeli,:);
    
            [~,pri2keep,~] = find(Cnull(vertcat(priCell{:}),:));
            pri2keep = unique(pri2keep);
            Cnull = Cnull(vertcat(priCell{:}),pri2keep);
    
            coarseSize(k,p,s,:) = numel(pri2keep);

            %% Compute initial condition from mp-problem
            [~, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load([nrbcell{:}]);
            msh_mp = msh_multipatch(mshCell,boundaries);
            spCurl_mp = sp_multipatch(spCurlCell,msh_mp,interfaces,boundary_interfaces);
            fInit = op_f_v_mp(spCurl_mp, msh_mp, initFun);
            ML2_mp = op_u_v_mp(spCurl_mp, spCurl_mp, msh_mp, @(x,y,z) ones(size(x)));
            [gDir_mp,dir_mp] = sp_drchlt_l2_proj(spCurl_mp, msh_mp, dirFunS, 1:numel(boundaries));
            vol_mp = setdiff(1:spCurl_mp.ndof,dir_mp);
    
            u0_mp = zeros(spCurl_mp.ndof,1);
            u0_mp(dir_mp) = gDir_mp*timeDep(t0);
            u0_mp(vol_mp) = ML2_mp(vol_mp,vol_mp)\(fInit(vol_mp) - ML2_mp(vol_mp,dir_mp)*u0_mp(dir_mp));
    
            u0Cell = cell(numel(nrbcell),1);
            for i=1:numel(nrbcell)
                u0Cell{i} = u0_mp(gnumCurl{i});
            end

            %% Minor preparations for class structure
            nPriCumuCell = cellfun(@(set) numel(set.pri), subsets, 'UniformOutput', false);
            nPriCumu = cumsum([0,nPriCumuCell{:}]);

            sys = FullSystemDD(subsets);
            sys.setCouplingRHS(zeros(numel(multNeli),1));
            sys.setNullRepresentations(cellfun(@(low,upp) Cnull(low:upp,:), num2cell(nPriCumu(1:end-1)+1), num2cell(nPriCumu(2:end)), 'UniformOutput', false));

            KCell = cellfun(@(sp, msh, ind) op_curlu_curlv_tp(sp, sp, msh, @(x,y,z) nu(x,y,z,ind)), spCurlCell, mshCell, indCell, 'UniformOutput', false);
            MCell = cellfun(@(sp, msh, ind) op_u_v_tp(sp, sp, msh, @(x,y,z) sigma(x,y,z,ind)), spCurlCell, mshCell, indCell, 'UniformOutput', false);
            fSCell = cellfun(@(sp, msh, ind) op_f_v_tp(sp, msh, @(x,y,z) sourceFunS(x,y,z,ind)), spCurlCell, mshCell, indCell, 'UniformOutput', false);

            for q = 1:numel(steps)
                fprintf('Simulating setup with N=%i, p=%i, s=%i, nT=%i:\n', numel(nrbcell),degrees(p),divs(s),steps(q));
                
                sys.reset();

                nT = steps(q);
                time = linspace(t0,tend,nT);
                dt = time(2) - time(1);

                %% Setup subdomain structure for time stepsize
                sys.setSystemMatrices(cellfun(@(M,K) M + dt*K, MCell, KCell, 'UniformOutput', false));
                sys.setCouplingMatrices(cellfun(@(low,upp) dt*B(:,low:upp), num2cell(ndofCumu(1:end-1)+1), num2cell(ndofCumu(2:end)), 'UniformOutput', false));
                sys.precomputeScaleMat(0);

                %% Compute coarse problem matrix
                sys.assembleCoarseF();
                sys.factorizeCoarseF();

                %% Evaluate RHS and eliminated function for every time step
                fTSCell = cellfun(@(f) dt*timeDep(time(2:end)).*f, fSCell, 'UniformOutput', false);
                geTSCell = cellfun(@(ge) timeDep(time(2:end)).*ge, geSCell, 'UniformOutput', false);

                %% Initialize solution
                solCell = u0Cell;

                %% Error measurement in every time step (with initial cond.)
                errBL2 = NaN*zeros(1,numel(time));
                errACL2 = NaN*zeros(1,numel(time));
                errECL2 = NaN*zeros(1,numel(time)-1);

                %% Measurements of other stuff
                condEstT = NaN*zeros(1,numel(time)-1);
                numIterT = NaN*zeros(1,numel(time)-1);

                for i = 1:numel(time)-1

                    %% Set RHS and eliminated term for time step
                    fCell = cellfun(@(fTS,M,x) fTS(:,i) + M*x(:,end), fTSCell, MCell, solCell, 'UniformOutput', false);
                    sys.setRHSs(fCell);
                    geCell = cellfun(@(geTS) geTS(:,i), geTSCell, 'UniformOutput', false);
                    sys.setEliminatedRHSs(geCell);

                    %% Computation of q
                    sys.computeOrthPriSol();

                    %% Computation of lamRem
                    sys.computeE();
                    sys.computeD();
                    [eigEst,numIterT(i)] = sys.computeLamRem(1e-6,1e3);
                    condEstT(i) = eigEst(2)/eigEst(1);

                    %% tRecovery of p
                    sys.recoverNullPriSol();

                    %% Recovery of ar
                    sys.recoverFullSol();
                    aCell = sys.getSol();

                    %% Compute error A and B
                    [~,errACell,errBCell] = cellfun(@(sp,msh,a) sp_hcurl_error(sp, msh, a, @(x,y,z) AFun(x,y,z,time(i+1)), @(x,y,z) BFun(x,y,z,time(i+1))), spCurlCell, mshCell, aCell, 'UniformOutput', false);
                    errBL2(i) = sqrt(sum([errBCell{:}].^2));
                    errACL2(i) = sqrt(sum([errACell{copper_regions}].^2));

                    %% Save solution for next time step
                    solCell = cellfun(@(sol,a) [sol,a], solCell, aCell, 'UniformOutput', false);

                    %% Compute error E
                    dta = (solCell{copper_regions}(:,i+1) - solCell{copper_regions}(:,i))/dt;
                    e = -dta;
                    [~,errECL2(i),~] = sp_hcurl_error(spCurlCell{copper_regions}, mshCell{copper_regions}, e, @(x,y,z) EFun(x,y,z,time(i+1)), @(x,y,z) BFun(x,y,z,time(i+1)));
                end

                %% Save measurements
                errBL2max(k,p,s,q) = max(errBL2);
                errACL2max(k,p,s,q) = max(errACL2);
                errECL2max(k,p,s,q) = max(errECL2);

                errB_acevedo(k,p,s,q) = sqrt(errBL2max(k,p,s,q)^2 + dt*sum(errBL2.^2));
                errE_acevedo(k,p,s,q) = sqrt(dt*sum(errECL2.^2));

                condEst(k,p,s,q) = mean(condEstT);
                numIter(k,p,s,q) = mean(numIterT);

                %% Export Measurements
                % subs, deg, divs, steps, pri, cond, iter, errB, errA, errE
                fprintf(fid, '%i,%i,%i,%i,%i,%.4d,%f,%.4d,%.4d,%.4d,%.4d,%.4d\n', numel(nrbcell), degrees(p), ...
                divs(s), nT, coarseSize(k,p,s,q), condEst(k,p,s,q), numIter(k,p,s,q), errBL2max(k,p,s,q),...
                errACL2max(k,p,s,q), errECL2max(k,p,s,q), errB_acevedo(k,p,s,q), errE_acevedo(k,p,s,q));

            end
        end
    end
end