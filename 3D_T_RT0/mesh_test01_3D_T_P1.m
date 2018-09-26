function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test01_3D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

elems2nodes = [1 2 3 4]';
nodes2coord = [0 0 0; 1 0 0; 0 1 0; 0 0 1]';

% maps for degrees of freedom (dof)
nelem = size(elems2nodes,2);
ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function 