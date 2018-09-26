function [K] = assembler2(Kloc,elems2faces,ndofs)
% function K = assemble(Kloc,elems2nodes,nodes2dofs,elist,ndof)
%Inefficient assembler, not getting uploaded
% ------------------------------------------------------------------------
% by Ted, March 2017
Agal = sparse(ndofs,ndofs);
sizes = size(elems2faces,2)*size(elems2faces,1)*size(elems2faces,1);
I = zeros(sizes,1);
J= zeros(sizes,1);
X = zeros(sizes,1);
ntriplets = 0;
for k = 1:size(elems2faces,2)
    for j = 1:size(elems2faces,1)
        for i = 1:size(elems2faces,1)
            %Agal(elems2faces(i,k),elems2faces(j,k)) = Agal(elems2faces(i,k),elems2faces(j,k)) + Kloc(i,j,k);
            ntriplets = ntriplets + 1;
            len = length (X) ;
            if (ntriplets > len)
                I (2*len) = 0 ;
                J (2*len) = 0 ;
                X (2*len) = 0 ;
            end
            I(ntriplets) = elems2faces(i,k);
            J(ntriplets) = elems2faces(j,k);
            X(ntriplets) = Kloc(i,j,k);
        end
    end
end
K = sparse(I,J,X,ndofs,ndofs);
%K = Agal;

return % end of function