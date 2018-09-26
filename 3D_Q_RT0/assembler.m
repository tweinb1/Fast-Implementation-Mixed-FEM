function [K] = assembler(Kloc,elems2faces,ndofs)
% function K = assemble(Kloc,elems2nodes,nodes2dofs,elist,ndof)
%Inefficient version of assembler
% Assemble a (global or subdomain) stiffness matrix. 
% Input: 
%   Kloc{e} ... local matrix of element (zeros at end if smaller),
%   elems2nodes{e}  ... global node number on element e,
%   nodes2dofs{a}   ... global dof numbers on global node a,
%   elist   ... list of elements to assemble; default=all.
%   ndof    ... (global) number of global dofs in the matrix K
% Output:
%   K       ... assembled stiffness matrix K.
% ------------------------------------------------------------------------
% by Ted, March 2017
Agal = sparse(ndofs,ndofs);
for k = 1:size(elems2faces,2)
    for j = 1:size(elems2faces,1)
        for i = 1:size(elems2faces,1)
            Agal(elems2faces(i,k),elems2faces(j,k)) = Agal(elems2faces(i,k),elems2faces(j,k)) + Kloc(i,j,k);
        end
    end
end
K = Agal;
%K = Agal;
system('free');
whos
return % end of function