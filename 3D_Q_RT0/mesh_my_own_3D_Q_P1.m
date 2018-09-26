function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_3D_Q_P1
% 
% Create mesh and setup degrees of freedom.
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% geometry of the physical domain
dim = 3;             % dimension of the physical space
dx = 2;              % dimension (size) of the domain in x-direction
dy = 2;              % dimension (size) of the domain in y-direction
dz = 2;
nx = 10; ny = 13; nz = 11;     % number of elements in one direction (on a side)
% --- end of set parameters --------------------------------------------

hx = dx/nx; hy = dy/ny;  hz = dz/nz; % mesh step 
% volD = dx*dy*dz;       % volume(=area) of the physical domain

[nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_3D_Q_P1(hx,hy,hz,nx,ny,nz);

% setup degrees of freedom
nodes2dofs = 1:(nx+1)*(ny+1); 
ngdof = max(max(nodes2dofs)); 

return % end of function