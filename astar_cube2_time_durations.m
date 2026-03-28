clear all;
close all;
restoredefaultpath();
addpath(genpath('./utils'))
addpath(genpath('./libraries'))

warning('off','geopdes:nrbmultipatch'); %-> Warns for using nrbmultipatch instead of mp_geo_load with files
warning('off','MATLAB:decomposition:LoadNotSupported'); %-> Composite structure don't allow to load decomposition objects (but are usable inside of spmd-blocks)

maxNumCompThreads(15); % Save resources if working on a Cluster

%% Setup numerical tests
deg = 3;
div = 2.^(3);
step = 2.^(5);
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

%% Setup of subdomains
nrbarr = cnstrct_cube_nrbarr([1,1,1],subs*[2,1,1]);
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


%% Setup for spaces and meshes
degree = deg*[1 1 1];
degCell = cellfun(@(nrb) degree, nrbcell, 'UniformOutput', false);
ndiv = div*[1 1 1];
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

nT = step;
time = linspace(t0,tend,nT);
dt = time(2) - time(1);

%% Setup measurement vectors
reps = 15;
durSetupCoarseMat = zeros(reps,1);
durSetupLocMat = zeros(reps,1);
durSetupRHS = zeros(reps,1);
durOrthPri = zeros(reps,1);
durCoarseRHS = zeros(reps,1);
durIntfPrb = zeros(reps,1);
durRecPri = zeros(reps,1);
durRecFull = zeros(reps,1);
durFull = zeros(reps,1);

for j=1:reps

    sys.reset();
    fprintf('Start run %i of %i\n',j,reps);

    tFull = tic;

    %% Setup subdomain structure for time stepsize
    tSetupLoc = tic;
    sys.setSystemMatrices(cellfun(@(M,K) M + dt*K, MCell, KCell, 'UniformOutput', false));
    sys.setCouplingMatrices(cellfun(@(low,upp) dt*B(:,low:upp), num2cell(ndofCumu(1:end-1)+1), num2cell(ndofCumu(2:end)), 'UniformOutput', false));
    sys.precomputeScaleMat(0);
    durSetupLocMat(j) = toc(tSetupLoc);
    
    %% Compute coarse problem matrix
    tSetupCoarseMat = tic;
    sys.assembleCoarseF();
    sys.factorizeCoarseF();
    durSetupCoarseMat(j) = toc(tSetupCoarseMat);
    
    %% Evaluate RHS and eliminated function for every time step
    tSetupRHS = tic;
    fTSCell = cellfun(@(f) dt*timeDep(time(2:end)).*f, fSCell, 'UniformOutput', false);
    geTSCell = cellfun(@(ge) timeDep(time(2:end)).*ge, geSCell, 'UniformOutput', false);
    durSetupRHS(j) = toc(tSetupRHS);
    
    %% Initialize solution
    solCell = u0Cell;
    
    for i = 1:numel(time)-1
    
        %% Set RHS and eliminated term for time step
        tSetupRHS = tic;
        fCell = cellfun(@(fTS,M,x) fTS(:,i) + M*x(:,end), fTSCell, MCell, solCell, 'UniformOutput', false);
        sys.setRHSs(fCell);
        geCell = cellfun(@(geTS) geTS(:,i), geTSCell, 'UniformOutput', false);
        sys.setEliminatedRHSs(geCell);
        durSetupRHS(j) = durSetupRHS(j) + toc(tSetupRHS);
    
        %% Computation of q
        tOrthPri = tic;
        sys.computeOrthPriSol();
        durOrthPri(j) = durOrthPri(j) + toc(tOrthPri);
    
        %% Computation of lamRem
        tCoarseRHS = tic;
        sys.computeE();
        sys.computeD();
        durCoarseRHS(j) = durCoarseRHS(j) + toc(tCoarseRHS);
    
        %% Solve interface problem
        tIntrfPrb = tic;
        sys.computeLamRem(1e-6,1e3);
        durIntfPrb(j) = durIntfPrb(j) + toc(tIntrfPrb);
    
        %% Recovery of p
        tRecPri = tic;
        sys.recoverNullPriSol();
        durRecPri(j) = durRecPri(j) + toc(tRecPri);
    
        %% Recovery of a
        tRecFull = tic;
        sys.recoverFullSol();
        durRecFull(j) = durRecFull(j) + toc(tRecFull);
    end
    durFull(j) = toc(tFull);
end

%% Plotting
figure(1)
clf()
timingsReal = [mean(durSetupLocMat), mean(durSetupCoarseMat), mean(durSetupRHS),...
    mean(durOrthPri), mean(durCoarseRHS), mean(durIntfPrb), mean(durRecPri),...
    mean(durRecFull)];
% Estimated performance in parallel implementation
timingsOpti = [mean(durSetupLocMat)/2, mean(durSetupCoarseMat), mean(durSetupRHS)/2,...
    mean(durOrthPri), mean(durCoarseRHS), mean(durIntfPrb), mean(durRecPri),...
    mean(durRecFull)/2];
bar({'Loc. Mats','Coarse Mat','Loc. RHSs','Orth. p','Coarse RHS','Intf. Prb.','Rec. p','Rec. a'},[timingsReal;timingsOpti]);
set(gca,'YScale','log');
grid on;
xtickangle(-65);
ylim([1e-3,1e1])
legend('Measured','Optimized')