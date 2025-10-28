function G = addEdgePrio2Graph_mortarIETI(G,dirSides,constrainedSides,freeSides,ietiSides)
    % dirSides = sides of dir boundary
    % constrainedSides = Constrained sides of mortaring interfaces
    % freeSides = free sides of mortaring interfaces
    % ietiSides = IETI and Periodicity BCs

    EdgeTable = G.Edges;
    NodeTable = G.Nodes;

    %% Iterate over all edges and assign priority
    geoInfo = EdgeTable.geoInfo;
    prio = zeros(numel(geoInfo),1);

    for i = 1:numel(geoInfo)
        prio(i) = compPrioIETI_dp(geoInfo{i},dirSides,constrainedSides,freeSides,ietiSides);
    end

    EdgeTable.Weight = prio;
    G = graph(EdgeTable,NodeTable);
end

function prio = compPrioIETI_dp(geoInfo,dir,constrained,free,ieti)
    
    cstr = union(ieti,constrained); % Priority right after Dirichlet
    not_nmn = union(dir,cstr);
    nmn = setdiff(1:6,not_nmn); % Implicitly contains free sides
    
    if numel(geoInfo)>=3 % check if some error occured
        error('Edges are not supposed to belong to more than two boundary sides!');
    elseif isempty(geoInfo) % Check if edge is in local volume
        prio = 7;
    elseif numel(geoInfo)==1 % check if edge does not belong to wireframe
        if nnz(ismember(geoInfo,nmn))==1 %Check if Neumann boundary
            prio = 6;
        elseif nnz(ismember(geoInfo,dir))==1 % Check if Dirichlet boundary
            prio = 6;
        elseif nnz(ismember(geoInfo,cstr))==1 % Check if constr interface
            prio = 6;
        else
            error('Found an edge with numel(geoInfo)==1 which is neither Dir., cstr. or Nmn.!');
        end
    elseif numel(geoInfo)==2 %Check if edge belongs to wireframe
        if nnz(ismember(geoInfo,dir))==2
            prio = 5; % Edge DD
        elseif nnz(ismember(geoInfo,dir))==1 && nnz(ismember(geoInfo,nmn))==1
            prio = 2; % Edge DN
        elseif nnz(ismember(geoInfo,dir))==1 && nnz(ismember(geoInfo,cstr))==1
            prio = 1; % Edge DC
        elseif nnz(ismember(geoInfo,nmn))==2
            prio = 5; % Edge NN
        elseif nnz(ismember(geoInfo,nmn))==1 && nnz(ismember(geoInfo,cstr))==1
            prio = 3; % Edge NC
        elseif nnz(ismember(geoInfo,cstr))==2
            prio = 4; % Edge CC
        else
            error('Something strange happened!')
        end
    end
end