function [nelem,ngnodes,elems2nodes,nodes2coord] = create_mesh_3D_Q_P1(hx,hy,hz,nx,ny,nz)
% 
% [nelem,nnelv,ien,id,xyz]
%(element,...
%    nx,ny,nz,dx,dy,dz)
%function [nelem,nnelv,ien,id,xyz]=create_mesh_hexahedr3d(element,...
%    nx,ny,nz,dx,dy,dz)
% Create mesh with hexahedral elements in 3D.
% ----------------------------------------------------------------------
% by Bedrich Sousedik, November 2015.

% number of finite elements
nelem = nx*ny*nz;

% coordinates
% [X,Y,Z] = meshgrid(hx*(0:nx),hy*(0:ny),hz*(0:nz));
[X,Y,Z] = meshgrid(hx*(0:nx),hy*(0:ny),hz*(0:nz));
X = permute(X,[2 1 3]); Y = permute(Y,[2 1 3]); Z = permute(Z,[2 1 3]);
nodes2coord = [X(:),Y(:),Z(:)]'; %,zeros((nx+1)*(ny+1),1)]';

% % coordinates
% hx = dx/nx; hy = dy/ny; hz = dz/nz;
% [X,Y,Z] = meshgrid(hx/2*(0:2*nx),hy/2*(0:2*ny),hz/2*(0:2*nz));
% xyz = [Y(:) X(:) Z(:)]; % Y first is not a typo


nnelv = 8;                         % local
ngnodes = (nx+1)*(ny+1)*(nz+1);    % global

% global number of node 1 on every element
xstep = 1;
ystep = 1*(nx+1);
%ystep = 1* (nx);
zstep = 1*(nx+1)*(ny+1);
%zstep = 1* (nx*ny);
node1gx = 1:xstep:nx;        % step in x direction
%node1gy = 1:ystep:nx*ny;     % step in y direction
node1gy = 1:ystep:ystep*ny;     % step in y direction
%node1gz = 1:zstep:nx*ny*nz;  % step in z direction
node1gz = 1:zstep:zstep*nz;  % step in z direction


node1g = zeros(nx,ny,nz);
layer1 = zeros(nx,ny);
layer1(1,:) = node1gy;
for i=1:ny
   layer1(:,i) = node1gy(i) + node1gx - 1;
end
for i=1:nz
   node1g(:,:,i) = layer1 + node1gz(i) - 1;  
end

% ien on the first element
ien1 = 1 + [ 0 1 ystep+1 ystep   ...     % 1-4
       zstep+[0 1 1+ystep ystep] ];%...     % 5-8
% ien1 = 1 + [ 0 2 ystep+2 ystep   ...     % 1-4
%        zstep+[0 2 2+ystep ystep] ...     % 5-8
%     1 2+1/2*ystep ystep+1 1/2*ystep ...  % 9-12
%   1/2*zstep+[0 2 2+ystep ystep] ...     % 13-16
%   zstep+[1 2+1/2*ystep ystep+1 1/2*ystep] ... %17-20
%     1/2*(xstep+ystep) 1/2*(xstep+zstep) 1/2*(ystep+zstep)+xstep ...
%     1/2*(xstep+zstep)+ystep 1/2*(ystep+zstep) 1/2*(xstep+ystep)+zstep ...
%     1/2*(xstep+ystep+zstep)];     
 
node1g = unique(node1g);

%Attempt edge case fix
if size(node1g,2) > 1
    node1g = node1g';
end
elems2nodes = (repmat(node1g,1,nnelv) + repmat(ien1,nelem,1) - 1)';


% % number of global pressure nodes : -------------------------------------   
% id = [ (1:ngnodes)' (ngnodes+(1:ngnodes))' (2*ngnodes+(1:ngnodes))' ...
%         zeros(ngnodes,1) ];
% %       [ (3*ngnodes+(1:npnodes))'; zeros(ngnodes-npnodes,1) ] ];
% % pressure dofs in the last column
% pnodes = unique(elems2nodes(1:8,:)); 
% id(pnodes,4) = 3*ngnodes+(1:nnz(pnodes));   
%   
% % id using PMD's numbering 1 2 3; 4 5; 6 7 8; ...
% % (this is used in create_sub.m in fields dof_node and dof_dof) 
% id2 = zeros(size(id));   pointer = 0;
% for i = 1:size(id,1) 
%    ndof = nnz(id(i,:)); 
%    id2(i,1:ndof) = pointer + (1:ndof); 
%    pointer = pointer + ndof;
% end    
% id = id2;
% % -----------------------------------------------------------------------


return % end of function