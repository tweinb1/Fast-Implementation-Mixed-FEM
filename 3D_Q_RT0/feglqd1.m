function[point,weight]=feglqd1(ngl)
% determine the integration points and weighting coefficients
% of Gauss-Legendre quadrature for one-dimensional integration
%-------------------------------------------------------------
% by Kannanut Chamsri and Bedrich Sousedik, June 2010.

% initialization
point  = zeros(ngl,1);
weight = zeros(ngl,1);

% find corresponding integration points and weights
switch ngl
   case 1   % 1-point quadrature rule
      point(1)  = 0;
      weight(1) = 2;
   case 2   % 2-point quadrature rule
      point(1)  = - 1/sqrt(3); %-0.577350269189626;
      point(2)  = 1/sqrt(3);  %-point(1);
      weight(1) = 1;
      weight(2) = 1;
   case 3   % 3-point quadrature rule
      point(1)  = - sqrt(3/5); %-0.774596669241483;
      point(2)  = 0;
      point(3)  = sqrt(3/5);  %-point(1);
      weight(1) = 5/9; %0.555555555555556;
      weight(2) = 8/9; %0.888888888888889;
      weight(3) = 5/9; %weight(1);   
   case 4    % 4-point quadrature rule   
      point(1)  = -0.861136311594053;
      point(2)  = -0.339981043584856;
      point(3)  = -point(2);
      point(4)  = -point(1);
      weight(1) = 0.347854845137454;
      weight(2) = 0.652145154862546;
      weight(3) = weight(2);  
      weight(4) = weight(1);  
   case 5    % 5-point quadrature rule   
      point(1)  = -0.906179845938664;
      point(2)  = -0.538469310105683;
      point(3)  = 0.0;
      point(4)  = -point(2);
      point(5)  = -point(1);
      weight(1) = 0.236926885056189;
      weight(2) = 0.478628670499366;
      weight(3) = 0.568888888888889;  
      weight(4) = weight(2);  
      weight(5) = weight(1);
   otherwise
      error('Unknown quadrature rule.')     
end
return % end of function   