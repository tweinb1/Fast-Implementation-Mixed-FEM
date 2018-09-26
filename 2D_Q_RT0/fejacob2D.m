function J = fejacob2D(dhdr,dhds,xcoord,ycoord)
% function J = fejacob2D(dhdr,dhds,xcoord,ycoord)
% Determine the Jacobian for two-dimensional mapping.
%  Variable Description:
%     dhdr - derivative of shape functions w.r.t. natural coordinate r
%     dhds - derivative of shape functions w.r.t. natural coordinate s
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
% ----------------------------------------------------------------------
% by Bedrich Sousedik, November 2015.
 
J = [ dhdr*xcoord' dhdr*ycoord'
      dhds*xcoord' dhds*ycoord' ];
    
return % end of function    