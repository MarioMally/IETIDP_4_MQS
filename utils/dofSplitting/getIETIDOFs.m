function [ieti] = getIETIDOFs(sp,msh,subIdx,intrfIETI)
    ietiPatches = [[intrfIETI.patch1],[intrfIETI.patch2]];
    ietiSides = [[intrfIETI.side1],[intrfIETI.side2]];

    if sp.ncomp==3
        [~,ieti] = sp_drchlt_l2_proj(sp,msh,@(x,y,z,i) zeros([3,size(x)]),ietiSides(ietiPatches==subIdx));
    elseif sp.ncomp==1
        [~,ieti] = sp_drchlt_l2_proj(sp,msh,@(x,y,z,i) zeros(size(x)),ietiSides(ietiPatches==subIdx));
    end
end