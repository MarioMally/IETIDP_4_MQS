function [msh,spCurl,spGrad,spDiv] = setupMeshAndSpaces(geo,nsub,deg,reg,nquad)
    [knots, zeta] = kntrefine(geo.nurbs.knots, nsub-1, deg, reg);
    [knots_hcurl, degree_hcurl] = knt_derham(knots, deg, 'Hcurl');
    rule = msh_gauss_nodes (nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh = msh_cartesian (zeta, qn, qw, geo);
    scalar_spaces = cell (msh.ndim, 1);
    for idim = 1:msh.ndim
        scalar_spaces{idim} = sp_bspline(knots_hcurl{idim}, degree_hcurl{idim}, msh);
    end
    spCurl = sp_vector(scalar_spaces, msh, 'curl-preserving');
    spGrad = sp_bspline(knots, deg, msh);

    if nargout==4
        [knots_hdiv, degree_hdiv] = knt_derham(knots, deg, 'Hdiv');
        scalar_spaces = cell (msh.ndim, 1);
        for idim = 1:msh.ndim
            scalar_spaces{idim} = sp_bspline(knots_hdiv{idim}, degree_hdiv{idim}, msh);
        end
        spDiv = sp_vector(scalar_spaces, msh, 'div-preserving');
    end
end