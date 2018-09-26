function [dhdx,dhdy]=federiv2(dhdr,dhds,invJ)
% function [dhdx,dhdy]=federiv2(dhdr,dhds,invJ)
% Determine derivatives of 2D isoparametric shape functions w.r.t.
% physical coordinate system.
%-----------------------------------------------------------------
% by Bedrich Sousedik, November 2015.

dhdx = invJ(1,1)*dhdr + invJ(1,2)*dhds;
dhdy = invJ(2,1)*dhdr + invJ(2,2)*dhds;

return % end of function