function [elems2edges, edges2nodes]=get_edges3(elems2nodes)
%function: [element2edges, edge2nodes]=get_edges(elems2nodes)
%requires: deleterepeatedrows
%generates edges of (triangular) triangulation defined in elems2nodes
%elems2nodes is matrix, whose rows contain numbers of its element nodes 
%element2edges returns edges numbers of each triangular element
%edge2nodes returns two node numbers of each edge
%example in 2D: [element2edges, edge2nodes]=get_edges([1 2 3; 2 4 3])
%example in 3D: [element2edges, edge2nodes]=get_edges([1 2 3 4; 1 2 3 5; 1 2 4 6])

%2D case
if (size(elems2nodes,2)==4)
    %extracts sets of edges 
    edges1=elems2nodes(:,[4 1]);
    edges2=elems2nodes(:,[1 2]);
    edges3=elems2nodes(:,[2 3]);
    edges4 = elems2nodes(:,[3 4]);
    alledges = [edges1;edges2;edges3;edges4];
    [test1,test2] = deleterepeatedrows(alledges);
    
    verticaledges = [edges1;edges3];
    horizontaledges = [edges2;edges4];
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