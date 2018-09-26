function [MASS,STIFF] = mass_pv_matrix_RT0(elems,B_K,B_K_det,signs,coeffs)
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
% Modified for Quadrilaterals December 2016


dim = 2;
nelems   = size(elems,1);        % number of elements
B_K_detA = abs(B_K_det);
% signs=abs(signs); % TEST!
%p = 1/sqrt(3);
%ip = [-p -p;
%       p -p;
%      -p  p;
%       p  p];

%[ip,w,nip] = intquad(2,dim);     % integration points and weighs (order 2, dimension 2 or 3)
%w = [1;1;1;1];
%nip = 4;

%Coeffs to reciprocal
coeffs = 1./coeffs;

[ip,w,nip] = intquadrilateral(2);
% reference basis function curl values on integration points,
% and number of basis functions
[val,dval,nbasis] = basis_RT0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
%STIFF = zeros(nbasis,nbasis,nelems);

STIFF = zeros(nbasis,1,nelems);
MASS = zeros(nbasis,nbasis,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis
            MASS(m,k,:) = squeeze(MASS(m,k,:))' + ...
                          w(i) .* B_K_detA'.^(-1) .* ... 
                          sum( squeeze( astam(signs(:,m), ( amsv(B_K, val(i,:,m)) .* reshape(coeffs,size(coeffs,1),1,size(coeffs,2)) ) ) ) ...
                               .* ...
                               squeeze( astam(signs(:,k), amsv(B_K, val(i,:,k))) ) ...
                             );
    end
        STIFF(m,1,:) = squeeze(STIFF(m,1,:)) + ...
                           w(i) .* B_K_detA.^(1) .* ...
                           ( signs(:,m) .* dval(i,:,m) );
    end
end

% copy symmetric entries of the local matrices
MASS = copy_triu(MASS);
% 
% % export
% MASS_loc = MASS;
% 
% % test of assembly
% %i=[1 1 2 2 2 2 3 3];j=i;
% %m(:,:,1)=ones(2);m(:,:,2)=m(:,:,1);
% %M=sparse(i,j,m(:))
% %y-indexes
% Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    
% % x-indexes (exchanges row and cols and preserves the third/elems dim.) 
% X = permute(Y,[2 1 3]); 
% % mass matrix
% MASS = sparse(X(:),Y(:),MASS(:)); 
% 
% % export
% STIFF_loc = STIFF;
% 
% %commented for now
% Y = elems';% Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    % y-indexes
% Y = Y(:);%permute(Y,[2 1 3]);                                       % x-indexes
% % new y indeces % ToDo
% X = ones(nbasis,nelems)*diag(1:nelems);
% STIFF = sparse(X(:),Y(:),STIFF(:));                           % stiffness matrix

end