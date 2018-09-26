function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_Matlab_3D_T_P1
% 
% 
% ----------------------------------------------------------------------
% by Bedrich Sousedik, May 2016.

dim = 3;

model = createpde(1);
% importGeometry(model,'BracketTwoHoles.stl');
importGeometry(model,'Block.stl');
generateMesh(model,'GeometricOrder','linear');% default: quadratic'););
pdeplot3D(model)

nelem = size(model.Mesh.Elements,2);

% elems2nodes : 
elems2nodes = model.Mesh.Elements;
% Hughes/my numbering of nodes 
% http://www.mathworks.com/help/pde/ug/finite-element-basis-for-3-d.html
% ç = pom;%zeros(size(pom)); 
% elems2nodes(1,:) = pom(2,:);
% elems2nodes(2,:) = pom(3,:);
% elems2nodes(3,:) = pom(4,:);
% elems2nodes(4,:) = pom(1,:);
% elems2nodes(5,:) = pom(6,:);
% elems2nodes(6,:) = pom(10,:);
% elems2nodes(7,:) = pom(8,:);
% elems2nodes(8,:) = pom(5,:);
% elems2nodes(9,:) = pom(9,:);
% elems2nodes(10,:) = pom(7,:);
% nelem=size(elems2nodes,2);ngnodes=size(nodes2coord,2);

nodes2coord = model.Mesh.Nodes;

nnodv = nnz(unique(model.Mesh.Elements)); % nodes with velocity dofs
nnodp = nnz(unique(model.Mesh.Elements(1:4,:))); % pressures only at triangle vertices
nodes2dofs = [ reshape(1:dim*nnodv,nnodv,dim)'; ...
               nnodv*dim+(1:nnodp) zeros(1,nnodv-nnodp) ];
ngdof = dim*nnodv + nnodp;

return % end of function