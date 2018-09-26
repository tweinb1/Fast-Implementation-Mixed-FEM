function [elems2faces, faces2nodes]=get_faces2(elems2nodes)
%function: [element2faces, face2nodes]=get_faces(elems2nodes)
%requires: deleterepeatedrows
%generates faces of (tetrahedral 3D) triangulation defined in elems2nodes
%elems2nodes is matrix, whose rows contain numbers of its element nodes 
%face2edges returns faces numbers of each tetrahedral element
%face2nodes returns three node numbers of each face
%example: [element2faces, face2nodes]=...
%get_faces([1 3 4 5; 7 4 3 5; 5 7 6 4; 6 8 7 4; 2 1 4 5; 2 4 6 5])
% ----------------------------------------------------------------------
% modified by Bedrich Sousedik, June 2016.

%extracts sets of faces
% faces1=elems2nodes(:,[1 2 3]);
% faces2=elems2nodes(:,[1 2 4]);
% faces3=elems2nodes(:,[1 3 4]);
% faces4=elems2nodes(:,[2 3 4]);
faces1=elems2nodes(:,[1 2 3 4]);
faces2=elems2nodes(:,[2 3 4 5]);
faces3=elems2nodes(:,[3 4 5 6]);
faces4=elems2nodes(:,[4 5 6 1]);
faces5=elems2nodes(:,[5 6 1 2]);
faces6=elems2nodes(:,[6 1 2 3]);

%as sets of their nodes (vertices)
vertices=zeros(size(elems2nodes,1)*6,4);
vertices(1:6:end,:)=faces1;
vertices(2:6:end,:)=faces2;
vertices(3:6:end,:)=faces3;
vertices(4:6:end,:)=faces4;
vertices(5:6:end,:)=faces5;
vertices(6:6:end,:)=faces6;


%repeated sets of nodes (joint faces) are eliminated 
[faces2nodes,elems2faces] = deleterepeatedrows(vertices);
elems2faces = reshape(elems2faces,6,size(elems2nodes,1))';

return % end of function