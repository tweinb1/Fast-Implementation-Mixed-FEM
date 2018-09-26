function [elems2edges, edges2nodes]=get_edges2(elems2nodes)
%function: [element2edges, edge2nodes]=get_edges2(elems2nodes)
%requires: deleterepeatedrows
%generates edges of quadrilaterals defined in elems2nodes
%elems2nodes is matrix, whose rows contain numbers of its element nodes 
%element2edges returns edges numbers of each quadrilateral element
%edge2nodes returns two node numbers of each edge
%By Theodore Weinberg, based on get_edges for triangles from Anjam

%2D case
if (size(elems2nodes,2)==4)
    %extracts sets of edges 
    edges1=elems2nodes(:,[2 3]);
    edges2=elems2nodes(:,[3 4]);
    edges3=elems2nodes(:,[4 1]);
    edges4 = elems2nodes(:,[1 2]);
    %Key for future:  
    %Edges1 gets you the rightmost vertical edge
    %edges2 gets the higher horizontal
    %edges3 gets left horizontal
    %edges4 gets lower vertical

    %as sets of their nodes (vertices)
    vertices=zeros(size(elems2nodes,1)*4,2);
    vertices(1:4:end,:)=edges1;
    vertices(2:4:end,:)=edges2;
    vertices(3:4:end,:)=edges3;
    vertices(4:4:end,:) = edges4;

    %repeated sets of nodes (joint edges) are eliminated 
    [edges2nodes,elems2edges]=deleterepeatedrows(vertices);
    elems2edges=reshape(elems2edges,4,size(elems2nodes,1))';
end