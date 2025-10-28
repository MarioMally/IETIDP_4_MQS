function locTree = glo2locTree(gloTree, gDOFs)
    lDOFs = 1:numel(gDOFs);
    locTree = lDOFs(ismember(gDOFs,gloTree));
end