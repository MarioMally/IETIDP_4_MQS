function [gloGraph] = loc2glo_graph(graphCell,gnum_curl,gnum_grad,offset)
    
    vec = cellfun(@(x) max(x), gnum_grad, 'UniformOutput', false);
    numGloNodes = max([vec{:}]);
    vec = cellfun(@(x) max(x), gnum_curl, 'UniformOutput', false);
    numGloEdges = max([vec{:}]);

    gloNodes = [1:numGloNodes]';
    gloEdges = zeros(numGloEdges,1);
    gloEndNodes = zeros(numGloEdges,2);
    gloWeight = inf*ones(numGloEdges,1);
    
    xCoor = inf*ones(numGloNodes,1);
    yCoor = inf*ones(numGloNodes,1);
    zCoor = inf*ones(numGloNodes,1);
%     geoInfo = [{},{}];

    for iPatch=1:numel(graphCell)
        G = graphCell{iPatch};
        % Node stuff
        gNodes = gnum_grad{iPatch};
        xCoor(gNodes) = G.Nodes.xCoor;
        yCoor(gNodes) = G.Nodes.yCoor;
        zCoor(gNodes) = G.Nodes.zCoor;
        

        % Edge stuff
        gEdges = gnum_curl{iPatch};
        ids = G.Edges.IDs;
        endNodes = str2double(G.Edges.EndNodes);

        gloEndNodes(gEdges(ids),:) = gNodes(endNodes);
        gloEdges(gEdges(ids)) =  gEdges(ids);

        for j=1:numel(ids)
            id = gEdges(ids(j));
            % Check for parabolic problem to identify intersection between conducting and nonconducting regions
            if (gloWeight(id)==2 &&  G.Edges.Weight(j)==offset + 2) || (gloWeight(id)==offset + 2 &&  G.Edges.Weight(j)==2)
                gloWeight(id) = 1;
            elseif (gloWeight(id)==4 &&  G.Edges.Weight(j)==offset + 4) || (gloWeight(id)==offset + 4 &&  G.Edges.Weight(j)==4)
                gloWeight(id) = 3;
            elseif (gloWeight(id)==7 &&  G.Edges.Weight(j)==offset + 7) || (gloWeight(id)==offset + 7 &&  G.Edges.Weight(j)==7)
                gloWeight(id) = 6;
            else
                gloWeight(id) = min([gloWeight(id),G.Edges.Weight(j)]);
            end
        end
    end
    
    %% Set up graph from tables
    NodeTable = table();
    EdgeTable = table();

    NodeTable.Name = cellfun(@(x) int2str(x), num2cell(gloNodes),'UniformOutput',false);
    NodeTable.IDs = gloNodes;
    NodeTable.xCoor = xCoor;
    NodeTable.yCoor = yCoor;
    NodeTable.zCoor = zCoor;
%     NodeTable.geoInfo = geoInfo';
%     NodeTable.Weight = gloNodePrio;

    EdgeTable.EndNodes = cellfun(@(x) int2str(x), num2cell(gloEndNodes),'UniformOutput',false);
    EdgeTable.IDs = gloEdges;
    EdgeTable.Weight = gloWeight;

    gloGraph = graph(EdgeTable,NodeTable);

end