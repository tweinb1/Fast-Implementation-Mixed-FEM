function K = assemble(Kloc,elems2faces,faces2dofs,elist,ndof)
% function K = assemble(Kloc,elems2nodes,nodes2dofs,elist,ndof)
% Assemble a (global or subdomain) stiffness matrix. 
% Input: 
%   Kloc{e} ... local matrix of element (zeros at end if smaller),
%   elems2faces{e}  ... global face number on element e,
%   faces2dofs{a}   ... global dof numbers on global face a,
%   elist   ... list of elements to assemble; default=all.
%   ndof    ... (global) number of global dofs in the matrix K
% Output:
%   K       ... assembled stiffness matrix K.
% ------------------------------------------------------------------------
% by Bedrich Sousedik, November 2015.

%tic

% preallocate
temp = cell(length(elist),1);

for e = elist(:)'
    
    % global node numbers or 0
    if iscell(elems2faces), 
       faces = elems2faces{e};
    else
       faces = elems2faces(:,e);
    end
    faces = faces(faces>0); % omit 0's
    
    % dofs on the global nodes or 0
    if iscell(faces2dofs), 
       dofs = [faces2dofs{faces}]; 
    else
       % dofs on the global nodes or 0
       dofs = faces2dofs(:,faces)';
    end
    dofs = dofs(dofs(:)>0); % stack grouped by node, omit 0's
    
    % take local as a 'sparse' matrix
    if iscell(Kloc), 
       [k,l,m] = find(Kloc{e});
    else
       [k,l,m] = find(Kloc(:,:,e));
    end
    % and insert it into the global and get indeces of the entries
    temp{e} = [dofs(k),dofs(l),m];
end

temp = cell2mat(temp);
K = sparse(temp(:,1),temp(:,2),temp(:,3),ndof,ndof);
%system('free');
%whos
%toc

return % end of function