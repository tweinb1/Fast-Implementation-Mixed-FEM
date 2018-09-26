function [Kloc,signs,JJ] = setup_stiffness_matrix_Darcy_2D_T_RT0(coeffs,nelem,elems2nodes,elems2edges,nodes2coord)
%function [Kloc,signs,J] = setup_stiffness_matrix_Darcy_2D_T_RT0(coeffs,nelem,elems2nodes,elems2edges,nodes2coord)
% Laplace equation on a triangle.
% local element numbering (cf. Hughes p. 168, here different):
%
%     3
%     | \   
%     1 - 2
%
% global numbering: x-directon first, then y-direction as 'layers'.
% ----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

dim = 2;

% number of dofs per element 
nelf = 3; % = 3 for 2D triangle, 4 for 3D tetrahedron
% number of dofs per element 
nedof = nelf + 1; % + pressure  

%Coeffs to reciprocal
coeffs = 1./coeffs;

% % just for check
% signs = signs_edges(elems2nodes'); signs=signs';
% signs0 = signs;
signs = zeros(size(elems2edges));

% quadrature for numerical integration (2*2 Gauss-Legendre quadrature)
ngl = 2;

% sampling points and weights
[point2,weight2] = feglqd2D_T(ngl);
% [point2,weight2,nip] = intquad(2,dim); weight2=2*weight2;

% computation of stiffness matrices and vectors and their assembly
Kloc = zeros(nedof,nedof,nelem); % init. element stiffness matrices

% disp('Entering the loop over elements (get matrices) ...')
for e = 1:nelem % loop over elements
   % signs of edges in this triangle
   tmp = elems2nodes([2 3 1],e) - elems2nodes([3 1 2],e);
   signs(:,e) = tmp ./ abs(tmp);
   
   % extract x, y values of the nodes
   xcoord = nodes2coord(1,elems2nodes(:,e));
   ycoord = nodes2coord(2,elems2nodes(:,e));
     
   % compute the Jacobian matrix
   J = [ xcoord(2)-xcoord(1) ycoord(2)-ycoord(1); ...
         xcoord(3)-xcoord(1) ycoord(3)-ycoord(1) ]  ;
   JJ(:,:,e)=J; 
   %invJ = inv(J);
   detJ = det(J)/2; % Jacobian %
   
   % and the element 'stiffness' coefficient
%    matmtx = coeffs(:,e);
   
   A = zeros(nelf);
   B = zeros(1,nelf);

   % numerical integration (loop over quadrature points)
   for ip = 1:size(point2,1) %ngl 
      x = point2(ip,1); 
      y = point2(ip,2);
      wip = weight2(ip);
      if abs(1-x-y) > eps, error(''), end
      
      % values of shape functions and derivatives at quadrature point
%       [sqhape,dhdr,dhds] = feisotriang2D3n(x,y);
      [sqhape,divsqhape] = feisotriang2D3n_RT0(x,y);
      
      % derivatives w.r.t. physical coordinate
%       [dhdx,dhdy] = federiv2D(dhdr,dhds,invJ);
      
      % compute the contribution in the element matrix
%        sqhape = diag(signs(:,e)) * sqhape * diag(coeffs(:,e));
      sqhape1 = spdiags(signs(:,e),0,nelf,nelf) * ( 1/(2*detJ) * (sqhape*J) * spdiags(coeffs(:,e),0,dim,dim) );
      sqhape2 = spdiags(signs(:,e),0,nelf,nelf) * 1/(2*detJ) * (sqhape*J);
      A = A + (sqhape1*sqhape2')*wip*detJ;
      B = B + ( signs(:,e)' .* divsqhape )*wip*detJ;
    %Ted:  moved Kloc line down
     % Kloc(:,:,e) = [A B'; 
                  %   B 0];
   end
     Kloc(:,:,e) = [A B'; 
                     B 0];
end

% % check
% norm(signs-signs0)

return % end of function