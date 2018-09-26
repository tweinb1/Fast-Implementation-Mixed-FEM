function signs_edges_Q = signs_edges_Q(elems2nodes)

% SIGNS_EDGES
% Returns signs for quadrilaterals
% March 2017 -Teddy
dim = size(elems2nodes,2)-2;

if ( dim == 2 )
    %tmp = elems2nodes(:,[2 3 4 1]) - elems2nodes(:,[3 4 1 2]);
    tmp = elems2nodes(:,[2 3 4 1]) - elems2nodes(:,[3 4 1 2]);
    signs_edges_Q = tmp ./ abs(tmp);
    %signs_edges_Q = ones(size(elems2nodes,1),size(elems2nodes,2));
    %signs_edges_Q(:,2:3) = -1;
    
elseif (dim == 3)
    tmp = elems2nodes(:,[1 1 1 2 3 4]) - elems2nodes(:,[2 3 4 3 4 2]);
    signs_edges_Q = tmp ./ abs(tmp);
else
    error('The data is not understood.')
end

end