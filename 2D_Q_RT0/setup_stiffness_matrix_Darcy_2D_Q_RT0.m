function [Kloc,signs] = setup_stiffness_matrix_Darcy_2D_Q_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord)
%
% ----------------------------------------------------------------------
% by ...
dim = 2;

%Coeffs to reciprocal
coeffs = 1./coeffs;

% number of dofs per element 
nelf = 4; % = 4 for 2D quadrilateral
% number of dofs per element 
nedof = nelf + 1; % + pressure  

% % just for check
% signs = signs_edges(elems2nodes'); signs=signs';
% signs0 = signs;
%signs = ones(size(elems2faces));
signs = signs_edges_Q(elems2nodes')';

%this point on unsure
% quadrature for numerical integration (2*2 Gauss-Legendre quadrature)
nglx = 2;
ngly=2;

% sampling points and weights
%Using ngl twice right now
[point2,weight2] = feglqd2D(nglx,ngly);
% [point2,weight2,nip] = intquad(2,dim); weight2=2*weight2;

% computation of stiffness matrices and vectors and their assembly
Kloc = zeros(nedof,nedof,nelem); % init. element stiffness matrices

%Ted Note:
%This point on only works for Triangles
% disp('Entering the loop over elements (get matrices) ...')
for e = 1:nelem % loop over elements
   % signs of edges in this triangle
   %Worrying about signs later for this moment, should be simpler for quad
  % tmp = elems2nodes([2 3 1],e) - elems2nodes([3 1 2],e);
   %signs(:,e) = tmp ./ abs(tmp);
   
   % extract x, y values of the nodes
   %Should work code as copied
   xcoord = nodes2coord(1,elems2nodes(:,e));
   ycoord = nodes2coord(2,elems2nodes(:,e));
   A = zeros(nelf);
   B = zeros(1,nelf);  
  for intx = 1:nglx
     x = point2(intx,1);
     wtx = weight2(intx,1);
     for inty = 1:ngly
         y = point2(inty,2);
         wty= weight2(inty,2);
         [sqhape,divsqhape,dhdr,dhds] = feisoquad2D4n_RT0(x,y);
         %Where to continue from
            % compute the Jacobian matrix
         
         %J = fejacob2D(dhdr,dhds,xcoord,ycoord);
         J = [(xcoord(2)-xcoord(1))/2 0;
             0 (ycoord(4)-ycoord(1))/2];
         %invJ = inv(J);
         detJ = det(J); % Jacobian
         
         % derivatives w.r.t. physical coordinates
        % [dhdx,dhdy] = federiv2(dhdr,dhds,invJ);
%trying without the 2
        %sqhape1 = spdiags(signs(:,e),0,nelf,nelf) * ( 1/(detJ) * (sqhape*J)); %* spdiags(coeffs(:,e),0,dim,dim) );
        sqhape1 = spdiags(signs(:,e),0,nelf,nelf) * ( 1/(detJ) * (sqhape*J)) * spdiags(coeffs(:,e),0,dim,dim);
        sqhape2 = spdiags(signs(:,e),0,nelf,nelf) * 1/(detJ) * (sqhape*J);
        A = A + (sqhape1*sqhape2')*wtx*wty*detJ;
        B = B + ( signs(:,e)' .* divsqhape )*wtx*wty*detJ;
         % compute the contribution in the element matrix
         %Kloc(:,:,e) = Kloc(:,:,e) + ...
          %   (dhdx'*matmtx*dhdx+dhdy'*matmtx*dhdy)*wtx*wty*detJ;
         
%          Mloc(:,:,e) = Mloc(:,:,e) + (sqhape'*sqhape)*wtx*wty*detJ;
         
         % ToDo: compute the right-hand side
         %f(:,1,e) = ...
         
     end
      Kloc(:,:,e) = [A B'; 
                     B 0];
  end
end
   % and the element 'stiffness' coefficient
%    matmtx = coeffs(:,e);
   


 
% % check
% norm(signs-signs0)
system('free');
whos;
          
return % end of function          