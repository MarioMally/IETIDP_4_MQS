function locGraph = setupLocalGraph(spCurl,nrb,subIdx,dirInfo,intrfsIETI,intfsMortar)
    locGraph = get_system_graph(spCurl,0,0,nrb);
    
    if ~isempty(dirInfo)
        dirPatches = [dirInfo.patch];
        dirSides = [dirInfo.side];
    else
        dirPatches = [];
        dirSides = [];
    end
    
    if ~isempty(intrfsIETI)
        ietiPatches = [[intrfsIETI.patch1],[intrfsIETI.patch2]];
        ietiSides = [[intrfsIETI.side1],[intrfsIETI.side2]];
    else
        ietiPatches = [];
        ietiSides = [];
    end

    if ~isempty(intfsMortar)
        constrPatches = [intfsMortar.constrPatch];
        constrSides = [intfsMortar.constrSide];
        freePatches = [intfsMortar.freePatch];
        freeSides = [intfsMortar.freeSide];
    else
        constrPatches = [];
        constrSides = [];
        freePatches = [];
        freeSides = [];
    end
    
    locGraph = addEdgePrio2Graph_mortarIETI(locGraph,...
        dirSides(dirPatches==subIdx), constrSides(constrPatches==subIdx),...
        freeSides(freePatches==subIdx), ietiSides(ietiPatches==subIdx));
end