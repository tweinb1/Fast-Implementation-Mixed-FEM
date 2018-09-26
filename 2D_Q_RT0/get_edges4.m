function [elems2edges, edges2nodes]=get_edges4(elems2nodes,numvert,numhoriz)

finaledges = zeros(size(elems2nodes,1)*4,1);
k = 1;
vertical = 0;
for i = 1:numvert
    for j = 1:4:numhoriz*4
        finaledges(j+vertical,1) = k;
        k = k+1;
    end
    k = k+1;
    vertical = vertical + numhoriz*4;
end


k = 2;
vertical = 0;
for i = 1:numvert
    for j=3:4:numhoriz*4
        finaledges(j+vertical) = k;
        k = k+1;
    end
    k = k+1;
    vertical = vertical + numhoriz*4;
end

vertical = 0;
k = (numhoriz+1) * numvert+1;
for i = 1:numvert
    for j = 2:4:numhoriz*4
        finaledges(j+vertical,1) = k;
        k = k+1;
    end
    vertical = vertical + numhoriz*4;
end

vertical = 0;
k = (numhoriz+1) * numvert + numhoriz + 1;
for i = 1:numvert
    for j = 4:4:numhoriz*4
        finaledges(j+vertical,1) = k;
        k = k + 1;
    end
    vertical = vertical + numhoriz*4;
end
elems2edges = finaledges;
elems2edges=reshape(elems2edges,4,size(elems2nodes,1))';