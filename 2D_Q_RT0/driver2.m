function [time, elengdof] = driver2(nx,ny)
%
% This is a generic driver for quadrilaterals
% This one does not do the non-vectorized
% nx: number of elements in x direction, nx >= 1
% ny: number of elements in y direction, ny >= 1
% ----------------------------------------------------------------------
% by Bedrich Sousedik and Theodore Weinberg, June 2016.

%tic
%tic
dim = 2;           % physical space dimension (only ndim = 2 supported)

% (mean) value of the coefficient

%Create mesh and show our figure, setup input data, currently random coeffs
coeffs_mean = 1;
[nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test_Q_P1(nx,ny);
[elems2faces,edges2nodes] = get_edges2(elems2nodes');
elems2faces = elems2faces';
nfaces = max(max(elems2faces));
elengdof = nfaces + nelem;
figure(1); show_mesh2(elems2nodes',nodes2coord'); title('mesh');%...
%coeffs = rand(dim,nelem)*10;
coeffs = ones(dim,nelem);
faces2nodes = edges2nodes;


%Vectorized
      
[jac,detj] = getDeterminants(elems2nodes,nodes2coord,dim,nelem);
signs = signs_edges_Q(elems2nodes');
%time(1) = toc;
tic
%mass_matrix_RT0(elems2faces,B_K,B_K_det,signs,coeffs)
[MPV] = mpv_matrix_RT0(elems2faces',jac,detj,signs,coeffs);
time(1) = toc;
system('free');
whos
tic
elems2nodes2 = [elems2faces; nfaces+(1:nelem)];
MPV0 = mat2cell(MPV,size(MPV,1),size(MPV,2),ones(size(MPV,3),1));
MPVA = assemble(MPV0,elems2nodes2,1:elengdof,1:nelem,elengdof);
time(2) = toc;
%tic
%MPVA2 = assembler(MPV, elems2nodes2, elengdof);
%time(4) = toc;
boundary =  [1,2,3,5,7,11,10,12];
bc = [boundary;0 0 0 0 0 0 0 0];
return %end of function