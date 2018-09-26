function [ ip, w, nip ] = intquad( qr, dim )

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
    [ip,w] = inttri(qr);
elseif (dim == 3)
    [ip,w] = inttet(qr);
else
    error('INTQUAD: Integration quadratures only for the unit triangle in 2D or unit tetrahedron in 3D.')
end

nip = size(ip,1);