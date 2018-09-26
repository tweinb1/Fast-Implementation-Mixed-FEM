function [nelem,elems2nodes,nodes2coord,nodes2dofs,ngdof] = mesh_Matlab_2D_T_P1
%
%
% ----------------------------------------------------------------------
% by Bedrich Sousedik, May 2016.

dim = 2;

model = createpde(1);
geometryFromEdges(model,@lshapeg);
generateMesh(model);%,'GeometricOrder','quadratic'); % default is linear
pdeplot(model)

nelem = size(model.Mesh.Elements,2);
elems2nodes = model.Mesh.Elements;
nodes2coord = model.Mesh.Nodes;

nnod = nnz(unique(model.Mesh.Elements)); 
nodes2dofs = 1:nnod;
ngdof = nnod;

return % end of function