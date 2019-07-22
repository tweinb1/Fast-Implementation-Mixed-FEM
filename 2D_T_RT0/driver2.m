function [time, elengdof] = driver2(nx,ny)
%
% This is a generic driver for quadrilaterals
% nx: number of elements in x direction, nx >= 1
% ny: number of elements in y direction, ny >= 1
% ----------------------------------------------------------------------
% by Bedrich Sousedik and Theodore Weinberg, June 2016.

%tic
dim = 2;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient

%Create mesh and show our figure, setup input data, currently random coeffs
coeffs_mean = 1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_num_2D_T_P1(nx,ny);
% adjust: facewise data for RT0 element
% note: regard edges in 2d as faces
elems2faces = get_edges(elems2nodes'); 
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
ngdof = nfaces + nelem;
elengdof = nfaces + nelem;
%figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...
%coeffs = rand(dim,nelem)*10;
coeffs = ones(dim,nelem);
%faces2nodes = edges2nodes;
elems2faces2 = [elems2faces; nfaces+(1:nelem)];

%Vectorized
[B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
signs         = signs_edges(elems2nodes');
tic
MPV = mpv_matrix_RT0(elems2faces',B_K,B_K_det,signs,coeffs); 
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
time(1) = toc;
tic
MPVA = assemble(MPV0,elems2faces2,1:ngdof,1:nelem,ngdof);
time(2) = toc;


return %end of function