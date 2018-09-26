function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test02_3D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% for 3d
load mesh3D

pom = elems2nodes';
elems2nodes = zeros(size(pom));
elems2nodes(1,:) = pom(4,:);
elems2nodes(2,:) = pom(3,:);
elems2nodes(3,:) = pom(2,:);
elems2nodes(4,:) = pom(1,:);
nodes2coord = nodes2coord';

% maps for degrees of freedom (dof)
nelem = size(elems2nodes,2);
ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function 