function [point3,weight3] = feglqd3(nglx,ngly,nglz)
% determine the integration points and weighting coefficients 
% of Gauss-Legendre quadrature for three-dimensinal integration.
%----------------------------------------------------------------
% by Kannanut Chamsri and Bedrich Sousedik, June 2010.

% initialization
ngl     = max([nglx ngly nglz]);
point3  = zeros(ngl,3);
weight3 = zeros(ngl,3);

% find corresponding integration points and weights for x,y,z - axes
[ponglx,weightx] = feglqd1(nglx);
[pongly,weighty] = feglqd1(ngly); 
[ponglz,weightz] = feglqd1(nglz); 
                                 
point3(1:nglx,1)  = ponglx(1:nglx);        
weight3(1:nglx,1) = weightx(1:nglx);
point3(1:ngly,2)  = pongly(1:ngly);
weight3(1:ngly,2) = weighty(1:ngly);
point3(1:nglz,3)  = ponglz(1:nglz);
weight3(1:nglz,3) = weightz(1:nglz);

return % end of function