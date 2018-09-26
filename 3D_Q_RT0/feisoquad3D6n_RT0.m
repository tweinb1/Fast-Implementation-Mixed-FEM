function [N,divN] = feisoquad3D6n_RT0(r,s,t)
% function [N,divN] = feisotetrahedr3D4n_RT0(r,s,t)
% Finite Element ISOparametric TETRAHEDRon in 3D with 4 nodes for RT0 elements
% Values of linear shape functions and their divergence. 
%-----------------------------------------------------------------------
% by Bedrich Sousedik, June 2016.

% values of the shape functions

N = .5 * [(1+r) 0 0;
            0 (1+s) 0;
            0 0 (1+t);
            (-1+r) 0 0;
            0 (-1+s) 0;
            0 0 (-1+t)];
        
% divergence
divN = zeros(1,6);
divN = divN + .5;
%divN(4:6) = -1 * divN(4:6);

return   % end of function