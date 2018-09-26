function [N,divN] = feisotetrahedr3D4n_RT0(r,s,t)
% function [N,divN] = feisotetrahedr3D4n_RT0(r,s,t)
% Finite Element ISOparametric TETRAHEDRon in 3D with 4 nodes for RT0 elements
% Values of linear shape functions and their divergence. 
%-----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% values of the shape functions
% val  = zeros( M , 3 , 4 );
%     val(:,:,1) = [ p(:,1)     p(:,2)     p(:,3)   ]*sqrt(3);
%     val(:,:,2) = [ p(:,1)-1   p(:,2)     p(:,3)   ];    
%     val(:,:,3) = [ p(:,1)     p(:,2)-1   p(:,3)   ];
%     val(:,:,4) = [ p(:,1)     p(:,2)     p(:,3)-1 ];
N = [ [  r    s    t ]*sqrt(3);
        r-1   s    t;
         r   s-1   t;
         r    s   t-1 ];

% divergence
% dval = zeros( M , 1 , 4 );
divN = zeros(1,4);
divN = divN + 3;
divN(1) = divN(1)*sqrt(3);

return   % end of function