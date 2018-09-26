function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_my_own_3D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% dx = 1/2;              % dimension (size) of the domain in x-direction
% dy = 1/2;              % dimension (size) of the domain in y-direction 
nx=3; ny=nx; nz=nx;     % number of elements in one direction (on a side)
% --- end of set parameters --------------------------------------------

% hx = dx/nx; hy = dy/ny;   % mesh size 
% volD = dx*dy;       % volume(=area) of the physical domain

% create mesh and setup degrees of freedom
[T,V] = regular_tetrahedral_mesh(nx+1,ny+1,nz+1);
tetramesh(T,V); % plot the mesh

elems2nodes = T';
nodes2coord = V';

pom = elems2nodes;
elems2nodes(1,:) = pom(2,:);
elems2nodes(2,:) = pom(3,:);
elems2nodes(3,:) = pom(4,:);
elems2nodes(4,:) = pom(1,:);

% maps for degrees of freedom (dof)
nelem = size(elems2nodes,2);
ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function