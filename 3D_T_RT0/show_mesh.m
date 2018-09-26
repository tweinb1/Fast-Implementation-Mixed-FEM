function show_mesh(elems2nodes,nodes2coord)
if (size(elems2nodes,2)==3)
    X=reshape(nodes2coord(elems2nodes',1),size(elems2nodes,2),size(elems2nodes,1));
    Y=reshape(nodes2coord(elems2nodes',2),size(elems2nodes,2),size(elems2nodes,1));
    %hold on
    fill(X,Y,[0.3 0.3 0.9]);

    % X=reshape(coordinates(elements4',1),size(elements4,2),size(elements4,1));
    % Y=reshape(coordinates(elements4',2),size(elements4,2),size(elements4,1));
    % fill(X,Y,[0.3 0.3 0.9]);
    % hold off
else
    tetramesh(elems2nodes,nodes2coord,'FaceAlpha',1);camorbit(20,0);
end


