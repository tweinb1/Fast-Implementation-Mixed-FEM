function [val,dval,nbasis] = basis_RT0(p)

% BASIS_RT0
%   The basis functions for linear Raviart-Thomas element for the
%   space H(div). This function evaluates these basis functions at
%   given points p in the reference geometry (quadrilateral in 2D and
%   block in 3D).
%
% Basis functions in 2D (edges):
%
% N = .5*[.5*(1 + r) 0;
%      0 .5*(1+s);
%      .5*(1-r) 0;
%       0 .5*(1-s)];

% Basis functions in 3D (faces):
%
% N = .125 * [(1+r) 0 0;
%             0 (1+s) 0;
%             0 0 (1+t);
%             (1-r) 0 0;
%             0 (1-s) 0;
%             0 0 (1-t)];
% SYNTAX:  [val,dval,nbasis] = basis_RT0(p)
%
%   In the following M denotes the number of integration poins.
%
% IN:   p       M x 2 or 3         vector defining the points
%
% OUT:  val     M x 2 or 3 x 3 or 4  basis function values on points p
%       dval    M x   1    x 3 or 4  basis function div values on points p
%       nbasis                       the number of basis functions
%
%       val(j,k,i): i'th basis function
%                   k'th component
%                   j'th point
%

%Updated for Quadrilaterals and Blocks by Theodore Weinberg


dim = size(p,2);

%2d
if ( dim == 2 )
    
    nbasis = 4;

    % initialize val and dval tensor
    M = size(p,1);
    val  = zeros( M , 2 , 4 );
    dval = zeros( M , 1 , 4 );

    % calculate basis function values
    val(:,:,1) = [ .5*(p(:,1)+1)     zeros(4,1)   ];
    val(:,:,2) = [ zeros(4,1)   .5*(1+p(:,2))];
    val(:,:,3) = [ .5*(-1+p(:,1))    zeros(4,1) ];
    val(:,:,4) = [ zeros(4,1)   .5*(-1+p(:,2))];

    
    % calculate basis function divergence values
    dval(:,:,1:2) = .5;
    dval(:,:,3:4) = .5;


%3d
elseif ( dim == 3 )
    
    nbasis = 6;
    
    % initialize val and dval tensor
    M = size(p,1);
    val  = zeros( M , 3 , 6 );
    dval = zeros( M , 1 , 6 );

    % calculate basis function values

    val(:,:,1) = [ p(:,1)+1    zeros(8,1) zeros(8,1)]*.5;
    val(:,:,2) = [ zeros(8,1)   p(:,2)+1     zeros(8,1)   ]*.5;    
    val(:,:,3) = [ zeros(8,1) zeros(8,1)    p(:,3)+1   ]*.5;
    val(:,:,4) = [ -1+p(:,1)     zeros(8,2) ]*.5;
    val(:,:,5) = [ zeros(8,1)   -1+p(:,2)     zeros(8,1)   ]*.5;    
    val(:,:,6) = [ zeros(8,2)     -1+p(:,3)   ]*.5;

    % calculate basis function divergence values
    dval(:,:,1:3) = .5;
    dval(:,:,4:6) = .5;
else
    
    error('BASIS_RT0: Input data not understood')
    
end