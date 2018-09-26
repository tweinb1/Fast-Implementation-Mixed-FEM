function [jac2,detj2] = getDeterminants(elems2nodes,nodes2coord,dim,nelem)
% function [jac2,detj2] = getDeterminants(elems2nodes,nodes2coord,dim,nelem)
% This will get the determinants and determinants of the jacobian in a
% format appropriate for the mass matrix rt0
% returns: jac2 = jacobian matrixx
%          detj2 = determinant of jacobian
%-----------------------------------------------------------------------
% by Theodore Weinberg and Bedrich Sousedik, December 2016.

%COULD USE UPDATE
p = 1/sqrt(3); %Value for quadrature
nlb = 4; %number of local basis functions, it must be known!
coord = zeros(dim,nlb,nelem);
for d = 1:dim
    for i = 1:nlb
          coord(d,i,:) = nodes2coord(d,elems2nodes(i,:));
    end
end

%Gets integration points
ip = [-p -p;
       p -p;
      -p  p;
       p  p]';

[dphi, detj, jac] = phider(coord,ip,'Q1');
jac = reshape(jac,2,2,nelem*4);
jac2 = zeros(2,2,nelem);
for i = 1:nelem
    jac2(:,:,i) = jac(:,:,4*i);
end

detj = abs(detj);
detj = reshape(detj,nelem*4,1);
detj2 = zeros(nelem,1);
for i = 1:nelem
    detj2(i,1) = detj(4*i,1);
end

return