function [Kloc] = mpv_matrix_RT0(elems,B_K,B_K_det,signs,coeffs)
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

%Coeffs to reciprocal
coeffs = 1./coeffs;

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
[ip,w,nip] = intquadrilateral(2);
% reference basis function curl values on integration points,
% and number of basis functions
[val,dval,nbasis] = basis_RT0(ip);

% calculate all local stiffness matrices simultaneously
% (without calculating symmetric entries twice)
%STIFF = zeros(nbasis,nbasis,nelems);

Kloc = zeros(nbasis+1,nbasis+1,nelems);
for i=1:nip
    for m=1:nbasis
        for k=m:nbasis
            Kloc(m,k,:) = squeeze(Kloc(m,k,:))' + ...
                          w(i) .* B_K_detA'.^(-1) .* ... 
                          sum( squeeze( astam(signs(:,m), ( amsv(B_K, val(i,:,m)) .* reshape(coeffs,size(coeffs,1),1,size(coeffs,2)) ) ) ) ...
                               .* ...
                               squeeze( astam(signs(:,k), amsv(B_K, val(i,:,k))) ) ...
                             );
         end
        Kloc(m,nbasis+1,:) = squeeze(Kloc(m,nbasis+1,:)) + ...
                           w(i) .* B_K_detA.^(1) .* ...
                           ( signs(:,m) .* dval(i,:,m) );
    end
end

% copy symmetric entries of the local matrices
Kloc = copy_triu(Kloc);
system('free');
whos
end