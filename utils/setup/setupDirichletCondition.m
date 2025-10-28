function [gDir,dirDOFs] = setupDirichletCondition(sp,msh,dirFun,subIdx,dirInfo)
    dirPatches = [dirInfo.patch];
    dirSides = [dirInfo.side];
    [gDir,dirDOFs] = sp_drchlt_l2_proj(sp,msh,dirFun,dirSides(dirPatches==subIdx));
end