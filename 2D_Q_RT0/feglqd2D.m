function [point2,weight2] = feglqd2D(nglx,ngly)
% function [point2,weight2] = feglqd2D(nglx,ngly)
% Determine the integration points and weighting coefficients 
% of Gauss-Legendre quadrature for two-dimensinal integration.
%-----------------------------------------------------------------------
% by Kannanut Chamsri and Bedrich Sousedik, June 2010.

% initialization
ngl     = max([nglx ngly]);
point2  = zeros(ngl,2);
weight2 = zeros(ngl,2);

% find corresponding integration points and weights for x,y,z - axes
[ponglx,weightx] = feglqd1(nglx);
[pongly,weighty] = feglqd1(ngly); 
                                 
point2(1:nglx,1)  = ponglx(1:nglx);        
weight2(1:nglx,1) = weightx(1:nglx);
point2(1:ngly,2)  = pongly(1:ngly);
weight2(1:ngly,2) = weighty(1:ngly);

return % end of function