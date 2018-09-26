function [Kloc,signs] = setup_stiffness_matrix_Darcy_3D_Q_RT0(...
             coeffs,nelem,elems2nodes,elems2faces,faces2nodes,nodes2coord)
%
% ----------------------------------------------------------------------
% by ...
dim = 3;

%Coeffs to reciprocal
coeffs = 1./coeffs;

% number of dofs per element
nelf = 6; % = 3 for 2D triangle, 4 for 3D tetrahedron
% number of dofs per element 
nedof = nelf + 1; % + pressure

% just for check
% [B_K,~,B_K_det] = affine_transformations(nodes2coord',elems2nodes');
% signs = signs_faces(nodes2coord',elems2faces',faces2nodes',B_K);
% signs0 = signs';
% signs = zeros(size(elems2faces));

% quadrature for numerical integration (2*2 Gauss-Legendre quadrature)
ngl = 2;
signs = ones(size(elems2faces));
% sampling points and weights
[point3,weight3] = feglqd3(ngl,ngl,ngl);
nglx = ngl;
ngly = ngl;
nglz = ngl;
% pre-allocate for matrices and vectors
% Kloc = cell(1,nelem);
Kloc = zeros(nedof,nedof,nelem); % init. element stiffness matrices
%floc = zeros(nedof,1,nelem); % init. element rhs vectors

% 1. Numerical quadrature
% ----------------------------------------------------------------------
% disp('Entering the loop over elements (get matrices) ...')
% computation of element matrices and vectors and their assembly
for e = 1:nelem
   % coordinates of the element nodes
   xcoord = nodes2coord(1,elems2nodes(:,e));
   ycoord = nodes2coord(2,elems2nodes(:,e));
   zcoord = nodes2coord(3,elems2nodes(:,e));
   

   % and the element 'stiffness' coefficient
%    matmtx = coeffs(e);
 
   % get signs for faces 
   ef = elems2faces(:,e);
%    for f = 1:length(ef)
%       p1 = nodes2coord(1:3,faces2nodes(1,ef(f)))';     % first,
%       p2 = nodes2coord(1:3,faces2nodes(2,ef(f)))';     % second,
%       p3 = nodes2coord(1:3,faces2nodes(3,ef(f)))';     % and third points of faces
%       vec1 = p2-p1;                               % two vectors defining
%       vec2 = p3-p1;                               % the face plane
%       normal = cross(vec1,vec2,2);               % normals
%        
%       n1 = [ 1  1  1]';
%       n2 = [-1  0  0]';
%       n3 = [ 0 -1  0]';
%       n4 = [ 0  0 -1]';
%       nf_ref = [n1 n2 n3 n4];
%       nf = J\nf_ref(:,f);
%       
%       signs(f,e) = sign(normal*nf);
%    end
   
   A = zeros(nelf);
   B = zeros(1,nelf);

    for intx = 1:nglx
      x = point3(intx,1);
      wtx = weight3(intx,1);
      for inty = 1:ngly
         y = point3(inty,2);
         wty = weight3(inty,2);
         for intz = 1:nglz
            z = point3(intz,3);
            wtz = weight3(intz,3);
               % compute the Jacobian matrix
            J = [(xcoord(2)-xcoord(1))/2 0 0;
                  0 (ycoord(3)-ycoord(1))/2 0;
                  0 0 (zcoord(5)-zcoord(1))/2];
            %    invJ = inv(J);
            %    JJ(:,:,e) = J;
              %detJ = det(J)/6; % Jacobian (1/6 to adjust for volume of tetrah., see Hughes p. 174)
             detJ = det(J);
        
                    % values of shape functions and derivatives at quadrature point
            [sqhape,divsqhape] = feisoquad3D6n_RT0(x,y,z);
%             [N,dNdksi,dNdeta,dNdzeta] = feisohexahedr3D8n(ksi,eta,zeta);               
%             
%             J = fejacob3D(dNdksi,dNdeta,dNdzeta,xcoord,ycoord,zcoord);
%             invJ = inv(J);  
%             detJ = det(J);
               
            % compute blocks of an element matrix    
         
               % compute the contribution in the element matrix
            sqhape1 = spdiags(signs(:,e),0,nelf,nelf) * ( 1/(detJ) * (sqhape*J) * spdiags(coeffs(:,e),0,dim,dim) );
            sqhape2 = spdiags(signs(:,e),0,nelf,nelf) * 1/(detJ) * (sqhape*J);
            A = A + (sqhape1*sqhape2')*wtx*wty*wtz*detJ;
            B = B + ( signs(:,e)' .* divsqhape )*wtx*wty*wtz*detJ;
            Kloc(:,:,e) = [A B'; 
                         B 0];
            % ToDo later (if needed):
            % element contribution to the rhs vector
            %floc = floc + Nv'*detJ*wtx*wty*wtz;
         end
      end
   end

      

      
      % derivatives w.r.t. physical coordinate
%       [dNdx,dNdy,dNdz] = federiv3D(dNvdksi,dNvdeta,dNvdzeta,invJ);

   

       % ToDo later (if needed):
       % element contribution to the rhs vector
       %floc = floc + Nv'*detJ*wtx*wty*wtz;
end


% B_K2 = permute(B_K,[2 1 3]); % Anjam-2015-FMA uses astually transpose of the Jacobian(!!!?)
% err = zeros(nelem,1); err2 = zeros(nelem,1);
% for i = 1:nelem%size(K,3)
%    err(i) = norm(JJ(:,:,i)-B_K2(:,:,i));
%    err2(i) = norm(JJ(:,:,i)-B_K(:,:,i)');
% end
% disp(sum(err)), disp(sum(err2))

%norm(signs-signs0)

system('free');
whos        
return % end of function          