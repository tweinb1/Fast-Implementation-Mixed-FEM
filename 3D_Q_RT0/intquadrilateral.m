function [ ip, w, nip ] = intquadrilateral(dim)

% INTQUAD
%   Integration quadrature for unit triangle in 2D or unit tetrahedron
%   in 3D. These quadrature rules provide means to exactly integrate
%   polynomials of given order.
%
% SYNTAX: [ ip, w, nip ] = intquad( qr, dim )
%
% IN:   qr      quadrature order (polynomials up to
%               order qr are integrated exactly)
%       dim     dimension (unit triangle in 2D, and
%               unit tetrahedron in 3D)
%
% OUT:  ip      integration points
%       w       integration weighs
%       nip     number of integration points/weighs
%

if ( dim == 2 )
    p = 1/sqrt(3);
    ip = [-p -p;
       p -p;
      -p  p;
       p  p];
   w = [1;1;1;1];
elseif (dim == 3)
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
    w = [1;1;1;1;1;1;1;1];
else
    error('Only for quadrilaterals or blocks')
end

nip = size(ip,1);