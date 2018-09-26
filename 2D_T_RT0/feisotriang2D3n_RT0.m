function [N,divN] = feisotriang2D3n_RT0(r,s)
% function [N,divN] = feisotriang2D3n_RT0(r,s)
% Finite Element ISOparametric Triangle in 2D with 3 nodes for RT0 elements
% Values of linear shape functions and their divergence. 
%-----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% values of the shape functions
% val  = zeros( M , 2 , 3 );
%     val(:,:,1) = [ p(:,1)     p(:,2)   ]*sqrt(2);
%     val(:,:,2) = [ p(:,1)-1   p(:,2)   ];
%     val(:,:,3) = [ p(:,1)     p(:,2)-1 ];
N = [ [  r    s  ]*sqrt(2);
        r-1   s  ;
         r   s-1   ];
    
% divergence
% dval = zeros( M , 1 , 3 );
divN = zeros(1,3);
divN = divN + 2;
divN(1) = divN(1)*sqrt(2);

return   % end of function