clear all;
close all;

[nrbarrCond2D, nrbarrAir2D, nrbarrCond3D, nrbarrAir3D] = planar_coil_mario(0.5, 0.6, 3, 0.5, 10, 20);

cond = 5.96*1e7;

%% Simulate steady state conduction
degree     = [3 3 3];  % Degree of the splines
regularity = [2 2 2];  % Regularity of the splines
nsub       = [8 8 8];  % Number of subdivisions
nquad      = [3 3 3];  % Points for the Gaussian quadrature rule

dirPatch1 = 1;
dirSide1 = 6;
dirPatch2 = 2;
dirSide2 = 6;

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (nrbarrCond3D);
npatch = numel (geometry);

vec1 = find([boundaries.patches]==dirPatch1 & [boundaries.faces]==dirSide1);
vec2 = find([boundaries.patches]==dirPatch2 & [boundaries.faces]==dirSide2);
dirSides = [vec1,vec2];

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

% Define the refined mesh, with tensor product structure
  [knots{iptc}, zeta{iptc}] = ...
         kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);

% Compute the quadrature rule
  rule      = msh_gauss_nodes (nquad);
  [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
  msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));

% Evaluate the discrete space basis functions in the quadrature points
  sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
end

msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp

% Compute and nrbarr_condassemble the matrices 
stiff_mat = op_gradu_gradv_mp (space, space, msh, @(x,y,z) ones(size(x))*cond);
rhs = op_f_v_mp (space, msh, @(x,y,z) zeros(size(x)));

% Apply Dirichlet boundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, @(x,y,z,ref) 10*(ref==dirSides(1))*ones(size(x)), dirSides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

vtk_pts = {linspace(0, 1, 30), linspace(0, 1, 30), linspace(0, 1, 30)};
sp_to_vtk(u, space, geometry, vtk_pts, 'data/conductor', {'phi','J'}, {'value','gradient'});
