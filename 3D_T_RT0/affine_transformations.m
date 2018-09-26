function [B_K,b_K,B_K_det] = affine_transformations(nodes2coord,elems2nodes)

% AFFINE_TRANSFORMATIONS
%   Calculate the affine transformations for a 2D triangular mesh, or
%   a 3D tetraherdal mesh, using the unit triangle/tetrahedron.
%
% SYNTAX: [B_K,b_K,B_K_det] = affine_transformations(node2coord,elems2nodes)
%
% IN:   nodes2coord     nodes by coordinates
%       elems2nodes     elements by nodes
%
% OUT:  B_K        the affine mapping from the unit triangle/tetrahedron:
%       b_K        F_K(x) = B_K * x + b_K
%       B_K_det    the determinant of B_K
%

dim = size(nodes2coord,2);

if ( dim == 2 )
    
    % initialize
    nelems = size(elems2nodes,1);
    B_K = zeros(2,2,nelems);

    % points defining the triangles
    A = nodes2coord(elems2nodes(:,1),1:2);
    B = nodes2coord(elems2nodes(:,2),1:2);
    C = nodes2coord(elems2nodes(:,3),1:2);
    
    % vectors defining the triangle
    a = B - A;
    b = C - A;
%     a = A - C;
%     b = B - C;

    % the affine mapping F_K = B_K * x + b_K
    % (bs: carefull - B_K is in fact a transpose of the Jacobian!)
    B_K(:,1,:) = a';
    B_K(:,2,:) = b';
    b_K        = A;

    % determinant
    B_K_det = a(:,1).*b(:,2) - a(:,2).*b(:,1);
    
elseif ( dim == 3 )
    
    % initialize
    nelems = size(elems2nodes,1);
    B_K = zeros(3,3,nelems);

    % points defining the elements
    A = nodes2coord(elems2nodes(:,1),1:3);
    B = nodes2coord(elems2nodes(:,2),1:3);
    C = nodes2coord(elems2nodes(:,3),1:3);
    D = nodes2coord(elems2nodes(:,4),1:3);
    
    % vectors defining the tetras
    a = B - A;
    b = C - A;
    c = D - A;
%     a = A - D;
%     b = B - D;
%     c = C - D;

    % the affine mapping F_K = B_K * x + b_K
    % this is in Anjam-2015-FMA:
    % (bs: carefull - B_K is in fact a transpose of the Jacobian!)
    B_K(:,1,:) = a';
    B_K(:,2,:) = b';
    B_K(:,3,:) = c';
    % bs this is correct, but their vectorization takes care of it
%     B_K(1,:,:) = a';
%     B_K(2,:,:) = b';
%     B_K(3,:,:) = c';
    b_K = A;
    
    % determinant
    cp = [ b(:,2).*c(:,3) - b(:,3).*c(:,2) ...    % crossproduct b x c
           b(:,3).*c(:,1) - b(:,1).*c(:,3) ...
           b(:,1).*c(:,2) - b(:,2).*c(:,1) ];
    B_K_det = sum( a .* cp, 2 );
    
else
    
    error('The input data is not understood.')
    
end