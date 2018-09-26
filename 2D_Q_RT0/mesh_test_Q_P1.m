function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_ted_test_2D_Q_P1(nx,ny)
% 
% Create mesh and setup degrees of freedom.
% nx: number of elements in x direction
% ny: number of elements in y direction
% ----------------------------------------------------------------------
% by Teddy Weinberg, December 2016.

% geometry of the physical domain
% ndim = 2;           % physical space dimension (only ndim = 2 supported)
dx = 2;              % dimension (size) of the domain in x-direction
dy = 2;              % dimension (size) of the domain in y-direction

% --- end of set parameters --------------------------------------------

hx = dx/nx; hy = dy/ny; % mesh step 
% volD = dx*dy;       % volume(=area) of the physical domain

[nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);

% setup degrees of freedom
nodes2dofs = 1:(nx+1)*(ny+1); 
ngdof = max(max(nodes2dofs)); 

return % end of function