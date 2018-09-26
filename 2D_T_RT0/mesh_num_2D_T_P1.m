function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_num_2D_T_P1(nx,ny)
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, May 2016.

dx = 1;              % dimension (size) of the domain in x-direction
dy = 1;              % dimension (size) of the domain in y-direction 
%nx = 2; ny = nx;     % number of elements in one direction (on a side)
% --- end of set parameters --------------------------------------------

hx = dx/nx; hy = dy/ny;   % mesh size 
% nelem = nx*ny;     % number of finite elements
% % volD = dx*dy;    % volume(=area) of the physical domain

% [elems2nodes,nodes2coord] = create_mesh_T_2D(hx,hy,nx,ny);
[nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_T_2D_P1(hx,hy,nx,ny);


% maps for degrees of freedom (dof)
% nelem = size(elems2nodes,2);
% ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function