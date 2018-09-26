function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_test01_2D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

elems2nodes = [3 1 2;1 3 4;5 1 4;1 5 6;7 1 6;1 7 8]'; %element;
elems2nodes = [2 3 1;4 1 3;4 5 1;6 1 5;6 7 1;8 1 7]'; %element;

% nodes2coord = 0.5*[0 0; 1 0; 2 0; 2 1; 1 1; 0 1; 0 2; 1 2]';
% elems2nodes = [1 2 5; 5 6 1; 2 3 4; 2 4 5; 6 5 8; 6 8 7]';
% dirichlet   = [1 2];

nodes2coord = [0 0; 1 0; 1 1; 0 1]';
elems2nodes = [2 3 1; 4 1 3]';

% nodes2coord = [0 0; 1 0; 0 1];
% elems2nodes = [1 2 3];

% maps for degrees of freedom (dof)
nelem = size(elems2nodes,2);
ngnodes = max(max(elems2nodes));
ngdof = ngnodes;   % number of global dofs
nodes2dofs = (1:ngdof);

return % end of function 