function [STIFF,STIFF_loc] = pv_matrix_RT0( elems, B_K_det, signs )
% based on
% function STIFF = stiffness_matrix_RT0( elems, B_K_det, signs )
% STIFFNESS_MATRIX_RT0
%   Vectorized calculation of stiffness matrix for the linear
%   Raviart-Thomas element in 2D/3D.
%
% SYNTAX:  STIFF = stiffness_matrix_RT0( elems, B_K_det, signs )
%
% IN:   elems    elements by: edges in 2D / faces in 3D
%       B_K_det  affine map determinants
%       signs    RT basis function signs per element
%
% OUT:  STIFF    the stiffness matrix
%

dim      = size(elems,2)-1;      % the dimension of the problem
nelems   = size(elems,1);        % number of elements
B_K_detA = abs(B_K_det);

p = 1/sqrt(3);
ip = [-p -p p;
       p -p p;
      -p  p p;
       p  p p;
      -p -p -p;
       p -p -p;
      -p  p -p;
       p  p -p
       ];
   
  
%NOTE:  PLEASE ADJUST THIS TO intquad format!!!!

%[ip,w,nip] = intquad(2,dim);     % integration points and weighs (order 2, dimension 2 or 3)
w = [1;1;1;1;1;1;1;1];
nip = 8;
% reference basis function curl values on integration points,
% and number of basis functions
[val,dval,nbasis] = basis_RT0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
%STIFF = zeros(nbasis,nbasis,nelems);
STIFF = zeros(nbasis,1,nelems);
for i=1:nip
    for m=1:nbasis
        for k=1:1%m:nbasis
%             STIFF(m,k,:) = squeeze(STIFF(m,k,:)) + ...
%                            w(i) .* B_K_detA.^(-1) .* ...
%                            ( signs(:,m) .* dval(i,:,m) );% .* ...
%                            %( signs(:,k) .* dval(i,:,k) );
            STIFF(m,k,:) = squeeze(STIFF(m,k,:)) + ...
                           w(i) .* B_K_detA.^(1) .* ...
                           ( signs(:,m) .* dval(i,:,m) );% .* ...
                           %( signs(:,k) .* dval(i,:,k) );
        end
    end
end

% copy symmetric entries of the local matrices
%STIFF = copy_triu(STIFF);

% export
STIFF_loc = STIFF;

Y = elems';% Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
Y = Y(:);%permute(Y,[2 1 3]);                                       % x-indexes
% new y indeces % ToDo
X = ones(nbasis,nelems)*diag(1:nelems);
STIFF = sparse(X(:),Y(:),STIFF(:));                           % stiffness matrix

end