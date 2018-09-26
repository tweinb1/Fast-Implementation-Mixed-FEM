function [nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny)
% function [nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_Q_2D_P1(hx,hy,nx,ny)
% Creates uniform mesh in 2D consisting of quadrilateral elements.
%  Input:
% hx,hy,nx,ny ... mesh sizes and number of elements in x, y directions
%  Output:
% nelem ... number of elements
% ngnodes ... number of {global) nodes
% elems2nodes ... element-node map (association)
% nodes2coord ... for each global node its x,y,z(=0) coordinates.
% ----------------------------------------------------------------------
% Bedrich Sousedik, November 2015.

% coordinates
[X,Y] = meshgrid(hx*(0:nx),hy*(0:ny)); X=X';Y=Y';
nodes2coord = [X(:),Y(:)]'; %,zeros((nx+1)*(ny+1),1)]';

% global node numbers on the grid
nodes = reshape(1:(nx+1)*(ny+1),nx+1,ny+1);
% global number of node 4 on every element
node4g = nodes(1:nx,1:ny);
node4 = node4g(:)';
elems2nodes = [ node4; node4+1; node4+nx+2; node4+nx+1];

nelem = nx*ny;             % number of finite elements
ngnodes = (nx+1)*(ny+1);   % number of global nodes

return % end of function