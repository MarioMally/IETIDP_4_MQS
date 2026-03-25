
function plot_graph_and_marked(G,marked,dir,shift,edgeNum,nodeNum,LW,NodeSize,cm)
    
    %% Plot graph and highlight marked
    H = plot(G,'EdgeLabel',edgeNum,'NodeLabel',nodeNum,'MarkerSize',0.1,'EdgeColor','k','NodeColor','k','EdgeLabelColor','k','NodeLabelColor','b','EdgeFontSize',10);
    H.XData = G.Nodes.xCoor+shift(1);
    H.YData = G.Nodes.yCoor+shift(2);
    H.ZData = G.Nodes.zCoor+shift(3);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
    H.EdgeAlpha = 1;
    
    if ~iscell(marked)
        T = graph();
        T = addnode(T,G.Nodes);
        for i=1:numel(marked)
            T = addedge(T,G.Edges(G.Edges.IDs==marked(i)+min(G.Edges.IDs)-1,:));
        end
        highlight(H,T,'EdgeColor','red','LineWidth',3);
    else
        % cm = distinguishable_colors(15,{'w','k','r'});
        % cm = ["#005aa9", "#ec6400", "#9ac103", "#710c84", "#000000"];
        % cm = ["#009D81", "#A60084"];
        for j=1:numel(marked)
            T = graph();
            T = addnode(T,G.Nodes);
            for i=1:numel(marked{j})
                T = addedge(T,G.Edges(G.Edges.IDs==marked{j}(i)+min(G.Edges.IDs)-1,:));
            end
            
            highlight(H,T,'EdgeColor',cm(j),'LineWidth',LW);
        end
    end

    %% Redraw nodes only for better graphics
    Gred = rmedge(G,G.Edges.IDs);
    H = plot(Gred,'NodeLabel',nodeNum,'MarkerSize',NodeSize,'NodeColor','k');
    H.XData = Gred.Nodes.xCoor+shift(1);
    H.YData = Gred.Nodes.yCoor+shift(2);
    H.ZData = Gred.Nodes.zCoor+shift(3);

    %% Mark specified local boundary sides
    faces = double.empty(0,4);
    v = [G.Nodes.xCoor'+shift(1);G.Nodes.yCoor'+shift(2);G.Nodes.zCoor'+shift(3)];

    %% Plot dependent interface sides
    for i=1:numel(dir)
        boolSide = cellfun(@(geoInfo) any(ismember(geoInfo,dir(i))),G.Nodes.geoInfo,'UniformOutput',false);
        boolSide = cell2mat(boolSide);
        subG = subgraph(G,find(boolSide));
        faceCell = allcycles(subG,"MaxCycleLength",4,"MinCycleLength",4);

        for j=1:numel(faceCell)
            face = str2double(faceCell{j});
            % Necessary to transform node index if a NodeOffset was used
            for k=1:4
                face(k) = find(G.Nodes.IDs==face(k));
            end
            faces = [faces;face];
        end
    end
    patch('Faces',faces,'Vertices',v','FaceColor','black','FaceAlpha',0.1,'EdgeColor','none');
    
    axis off;
    %% Plot independent interface sides
%     for i=1:numel(idep)
%         boolSide = cellfun(@(geoInfo) any(ismember(geoInfo,idep(i))),G.Nodes.geoInfo,'UniformOutput',false);
%         boolSide = cell2mat(boolSide);
%         subG = subgraph(G,find(boolSide));
%         faceCell = allcycles(subG,"MaxCycleLength",4,"MinCycleLength",4);
% 
%         for j=1:numel(faceCell)
%             face = str2double(faceCell{j});
%             % Necessary to transform node index if a NodeOffset was used
%             for k=1:4
%                 face(k) = find(G.Nodes.IDs==face(k));
%             end
%             faces = [faces;face];
%         end
%     end
%     patch('Faces',faces,'Vertices',v','FaceColor','yellow','FaceAlpha',.2,'EdgeColor','none');
    
end