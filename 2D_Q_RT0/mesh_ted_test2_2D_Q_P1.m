function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_ted_test_2D_Q_P1
% 
% Create mesh and setup degrees of freedom.
% ----------------------------------------------------------------------
% by Teddy Weinberg, June 2016.

% geometry of the physical domain
% ndim = 2;           % physical space dimension (only ndim = 2 supported)
dx = 1;              % dimension (size) of the domain in x-direction
dy = 2;              % dimension (size) of the domain in y-direction
nx = 3; ny = 3;%nx;     % number of elements in one direction (on a side)
% --- end of set parameters --------------------------------------------

hx = dx/nx; hy = dy/ny; % mesh step 
% volD = dx*dy;       % volume(=area) of the physical domain

[nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny);

% setup degrees of freedom
nodes2dofs = 1:(nx+1)*(ny+1); 
ngdof = max(max(nodes2dofs)); 

return % end of function