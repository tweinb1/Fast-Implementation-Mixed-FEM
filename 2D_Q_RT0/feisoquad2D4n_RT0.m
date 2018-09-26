function [N,divN,dhdrs,dhdss] = feisoquad2D4n_RT0(r,s)
% function [N,divN] = feisoquad2D3n_RT0(r,s)
% Based on Finite Element ISOparametric Quadrilateral in 2D with 4 nodes for RT0 elements
% Based on Triangle code
% Values of linear shape functions and their divergence. 
%-----------------------------------------------------------------------
% by Theodore Weinberg and Bedrich Sousedik, June 2016.

% values of the shape functions
% val  = zeros( M , 2 , 3 );
%     val(:,:,1) = [ p(:,1)     p(:,2)   ]*sqrt(2);
%     val(:,:,2) = [ p(:,1)-1   p(:,2)   ];
%     val(:,:,3) = [ p(:,1)     p(:,2)-1 ];
%N = [ [  r    s  ]*sqrt(2);
     %   r-1   s  ;
       %  r   s-1   ];
N = [.5*(1 + r) 0;
     0 .5*(1+s);
     .5*(-1+r) 0;
      0 .5*(-1+s)];

 dhdrs = [.5 0;
          0 0;
          .5 0;
          0 0]';
 dhdss = [0 0;
          0 .5;
          0 0;
          0 .5]';
% divergence
% dval = zeros( M , 1 , 3 );
divN = zeros(1,4);
divN(1) = .5;
divN(2) = .5;
divN(3) = .5;
divN(4) = .5;
%divN(1) = divN(1)*sqrt(2);

return   % end of function