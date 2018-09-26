function [point2,weight2] = feglqd2D_T(ngl)
% function [point2,weight2] = feglqd2D_T(ngl)
% Determine the integration points and weighting coefficients of 
% Gauss-Legendre quadrature for two-dimensinal integration on a triangle.
% Used Table 3.I.1, p. 173 of Hughes FE-book. 
% ------------------------------------------------------------------------
% by Bedrich Sousedik, November 2015.

% find corresponding integration points and weights
switch ngl
    case 1
        point2 = [1/3 1/3];
        weight2 = 1;
    case 2 % 3-point formula, degree of precision 2
        point2 =  [ 1/2  0 ;% 1/2 ;   % point 1
                     0  1/2;% 1/2;   % point 2
                    1/2 1/2];%  0 ];  % point 3
        weight2 = [ 1/3 1/3 1/3];
    case 3 % 4-point formula, % degree of precision 3
        point2 = [ 1/3 1/3;
                   3/5 1/5;
                   1/5 1/5;
                   1/5 3/5];
        weight2 = [-27/48 25/48*ones(1,3)];      
%         point2 = [ 1/3 1/3 1/3;   % point 1
%                    0.6 0.2 0.2;   % point 2
%                    0.2 0.6 0.2;   % point 3
%                    0.2 0.2 0.6 ]; % point 4
%         weight2 = [ -0.56250 0.520833333333333 * ones(1,3) ];
%    case 3 % 6-point formula, degree of precision 3
%         p1 = 0.659027622374092;
%         p2 = 0.231933368553031; 
%         p3 = 0.109039009072877;
%         point2 = [ p1 p2 p3;   % point 1
%                    p3 p1 p2;   % point 2
%                    p2 p3 p1;   % point 3
%                    p3 p2 p1;   % point 4
%                    p1 p3 p2;   % point 5
%                    p2 p1 p3];  % point 6
%         weight2 = 1/6 * ones(1,6);
    case 4 % 6-point formula, degree of precision 4
        p1 = 0.816847572980459;
        p2 = 0.091576213509771;
        p3 = 0.108103018168070;
        p4 = 0.445948490915965;
        w1 = 0.109951743655322;
        w2 = 0.223381589678011;
        point2 = [ p1 p2 p2;   % point 1
                   p2 p1 p2;   % point 2
                   p2 p2 p1;   % point 3
                   p3 p4 p4;   % point 4
                   p4 p3 p4;   % point 5
                   p4 p4 p3 ]; % point 6
        weight2 = [ w1 w1 w1 w2 w2 w2 ];   
    otherwise
      error('Unknown quadrature rule.')     
end

return % end of function