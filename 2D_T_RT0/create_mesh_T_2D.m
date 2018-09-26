function [elems2nodes,nodes2coord] = create_mesh_T_2D(hx,hy,nx,ny)
% function [elems2nodes,nodes2coord] = create_mesh_T_2D(hx,hy,nx,ny)
%  input:
% hx,hy,nx,ny ... mesh sizes and number of elements in x, y directions
%  output:
% elems2nodes ... element-node map (association)
% nodes2coord ... for each global node its x,y,z(=0) coordinates.
% ----------------------------------------------------------------------
% Bedrich Sousedik, November 2015.

% coordinates
% hx=1;hy=1/2;nx=4;ny=2; % for testing
[X,Y] = meshgrid(hx*(0:nx),hy*(0:ny)); X=X';Y=Y';
nodes2coord = [X(:),Y(:)]; %,zeros((nx+1)*(ny+1),1)];
elems2nodes = delaunay(X,Y);

return % end of function