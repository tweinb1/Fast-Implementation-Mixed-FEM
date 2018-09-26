function [MASS,MASS_loc] = mass_matrix_RT0(elems,B_K,B_K_det,signs,coeffs)

% MASS_MATRIX_RT0
%   Vectorized calculation of mass matrix for the linear
%   Raviart-Thomas element in 2D/3D.
%
% SYNTAX:  MASS = mass_matrix_RT0( elems, B_K, B_K_det, signs )
%
% IN:   elems    elements by: edges in 2D / faces in 3D
%       B_K      affine map matrices
%       B_K_det  affine map determinants
%       signs    RT0 basis function signs per element
%
% OUT:  MASS     the mass matrix
%

dim      = size(elems,2)-1;      % the dimension of the problem
nelems   = size(elems,1);        % number of elements
B_K_detA = abs(B_K_det);
% signs=abs(signs); % TEST!

[ip,w,nip] = intquad(2,dim);     % integration points and weighs (order 2, dimension 2/3)

% reference basis function curl values on integration points,
% and number of basis functions
[val,~,nbasis] = basis_RT0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
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
    end
end

% copy symmetric entries of the local matrices
MASS = copy_triu(MASS);

% export
MASS_loc = MASS;

% test of assembly
%i=[1 1 2 2 2 2 3 3];j=i;
%m(:,:,1)=ones(2);m(:,:,2)=m(:,:,1);
%M=sparse(i,j,m(:))

% y-indexes
Y = reshape(repmat(elems',nbasis,1),nbasis,nbasis,nelems);    
% x-indexes (exchanges row and cols and preserves the third/elems dim.) 
X = permute(Y,[2 1 3]); 
% mass matrix
MASS = sparse(X(:),Y(:),MASS(:)); 

% added by bs - using the loop over elements
% M2 = zeros(nbasis,nbasis,nelems);
% for e=1:size(elems,1)
%    % loop over quadrature points   
%  
%    M2(i,j,e)=1; 
% end

end