function [time, ngdof] = driver2 (nx, ny, nz)
%  This is driver, except it only does vectorized code
%
%
% ----------------------------------------------------------------
% by Theodore Weinberg, November 2016

dim = 3;

% (mean) value of the coefficient
coeffs_mean = 1;

[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_num_3D_T_P1(nx,ny,nz);

% adjust: facewise data for RT0 element
% ToDo: elems2faces = get_edges(elems2nodes)
% etc.
[elems2faces,faces2nodes] = get_faces(elems2nodes');
elems2faces = elems2faces'; 
faces2nodes = faces2nodes';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;
coeffs = 1*ones(dim,nelem);
coeffs = rand(dim,nelem);
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
%...
% affine transformations (same 2d and 3d, but only for triangles or tetrahedrons)
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
signs = signs_faces(nodes2coord',elems2faces',faces2nodes',B_K);

tic;
[MPV] = mpv_matrix_RT0(elems2faces',B_K,B_K_det,signs,coeffs);
time(1) = toc;
tic;
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2nodes2,1:ngdof,1:nelem,ngdof);
time(2) = toc;



return %end of function