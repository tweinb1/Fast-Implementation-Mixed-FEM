function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test02_2D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% for 2d
load testtest
elems2nodes = element';
nodes2coord = coordinate';

% maps for degrees of freedom (dof)
nelem = size(elems2nodes,2);
ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function 