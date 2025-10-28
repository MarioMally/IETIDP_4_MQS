function C = patch2multipatch(gnum)
    gDOFs = []; lDOFs = [];
    cumuNdofs = 0;
    for i=1:numel(gnum)
        ndof = numel(gnum{i});
        dofs = 1:ndof;
        gDOFs = [gDOFs,gnum{i}];
        lDOFs = [lDOFs,cumuNdofs + dofs];
        cumuNdofs = cumuNdofs + ndof;
    end
    gDOFs_unique = unique(gDOFs);

    glDOFs  = gDOFs(lDOFs);

    n = cumuNdofs;
    m = numel(gDOFs_unique);
    [~,rows] = ismember(glDOFs,gDOFs_unique);

    C = sparse(rows,lDOFs,ones(numel(glDOFs),1),m,n)';
    %Remove zero rows
    C = C(:,ismember(1:size(C,2),rows));
end