function [val,dval,nbasis] = basis_RT0(p)

% BASIS_RT0
%   The basis functions for linear Raviart-Thomas element for the
%   space H(div). This function evaluates these basis functions at
%   given points p in the reference geometry (triangle in 2D and
%   tetrahedron in 3D).
%
% Basis functions in 2D (edges):
%
%   N_1(x,y) = [ x ]     N_2(x,y) = [ x - 1 ]     N_3(x,y) = [   x   ]
%              [ y ]                [   y   ]                [ y - 1 ]
%
% Basis functions in 3D (faces):
%
%   N_1(x,y,z) = [   x   ]      N_2(x,y,z) = [   x   ]
%                [   y   ]                   [ y - 1 ]
%                [ z - 1 ]                   [   z   ]
%
%   N_3(x,y,z) = [ x - 1 ]      N_4(x,y,z) = [ x ]
%                [   y   ]                   [ y ]
%                [   z   ]                   [ z ]
%
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

dim = size(p,2);

%2d
if ( dim == 2 )
    
    nbasis = 3;

    % initialize val and dval tensor
    M = size(p,1);
    val  = zeros( M , 2 , 3 );
    dval = zeros( M , 1 , 3 );

    % calculate basis function values
    val(:,:,1) = [ p(:,1)     p(:,2)   ]*sqrt(2);
    val(:,:,2) = [ p(:,1)-1   p(:,2)   ];
    val(:,:,3) = [ p(:,1)     p(:,2)-1 ];
%     val(:,:,1) = [ p(:,1)-1   p(:,2)   ];
%     val(:,:,2) = [ p(:,1)     p(:,2)-1 ];
%     val(:,:,3) = [ p(:,1)     p(:,2)   ]*sqrt(2);
    
    % calculate basis function divergence values
    dval = dval + 2;
    dval(:,:,1) = dval(:,:,1)*sqrt(2); 
%     dval(:,:,3) = dval(:,:,3)*sqrt(2);

%3d
elseif ( dim == 3 )
    
    nbasis = 4;
    
    % initialize val and dval tensor
    M = size(p,1);
    val  = zeros( M , 3 , 4 );
    dval = zeros( M , 1 , 4 );

    % calculate basis function values
%     val(:,:,1) = [ p(:,1)     p(:,2)     p(:,3)-1 ];
%     val(:,:,2) = [ p(:,1)     p(:,2)-1   p(:,3)   ];
%     val(:,:,3) = [ p(:,1)-1   p(:,2)     p(:,3)   ];
%     val(:,:,4) = [ p(:,1)     p(:,2)     p(:,3)   ];
    val(:,:,1) = [ p(:,1)     p(:,2)     p(:,3)   ]*sqrt(3);
    val(:,:,2) = [ p(:,1)-1   p(:,2)     p(:,3)   ];    
    val(:,:,3) = [ p(:,1)     p(:,2)-1   p(:,3)   ];
    val(:,:,4) = [ p(:,1)     p(:,2)     p(:,3)-1 ];

    % calculate basis function divergence values
    dval = dval + 3;
    dval(:,:,1) = dval(:,:,1)*sqrt(3);
else
    
    error('BASIS_RT0: Input data not understood')
    
end