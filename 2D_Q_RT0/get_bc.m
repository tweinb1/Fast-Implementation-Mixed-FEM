function bc = get_bc(edges2dofs,boundary)
%GET_BC_3DHASAL
% Get boundary conditions for 3D T Hasal experiment
%
% Input:
%   nodes2dofs  - mapping of nodes to degrees of freedom
%   nodes2coord - mapping of nodes to coordinates
%   boundary    - boundary information from netgen
%
% Output:
%   bc - boundary conditions

% compile boundary nodes
%bnodes = boundary(2:7,:);
%bnodes = unique(bnodes(:));
%nnodes = nnz(bnodes);
nnodes = nnz(boundary);

% initialize dof vector and value vector
bcdof  = edges2dofs(1:4,boundary);
bcval  = zeros(3,nnodes);


% loop over nodes and set boundary condition
for node = 1:nnodes
    bcval(2,node) = 0;
end

% get rid of Neumann do-nothing nodes
bcdof = bcdof(:)';
bcval = bcval(:)';

% reorder (to have dofs in increasing/ascending order)
[bcdof,index] = sort(bcdof,'ascend');
bcval = bcval(index); 
bc = [bcdof; bcval];
end