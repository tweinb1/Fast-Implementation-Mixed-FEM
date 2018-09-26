function [K] = assembler(Kloc,elems2faces,ndofs)
% function K = assemble(Kloc,elems2nodes,nodes2dofs,elist,ndof)
%Inefficient version of assembler
% Assemble a (global or subdomain) stiffness matrix. 
% Input: 
%   Kloc ... local matrices
%   elems2faces  ... corresponding elements to their edges/faces,
%   ndofs    ... (global) number of global dofs in the matrix K
% Output:
%   K       ... assembled stiffness matrix K.
% ------------------------------------------------------------------------
% by Ted, March 2017
K = sparse(ndofs,ndofs);
for k = 1:size(elems2faces,2)
    for j = 1:size(elems2faces,1)
        for i = 1:size(elems2faces,1)
            K(elems2faces(i,k),elems2faces(j,k)) = K(elems2faces(i,k),elems2faces(j,k)) + Kloc(i,j,k);
        end
    end
end
%system('free');
%whos
return % end of function