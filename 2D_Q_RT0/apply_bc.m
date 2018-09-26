function [A,f] = apply_bc(A,f,bc)
% 
% Apply boundary conditions to global matrix A and the rhs vector f.
%  Input:
% A ... system matrix 
% f ... system rhs
% bcdof ... a node vector of constrained dofs
% bcval ... for each constrained dof its value.
% ----------------------------------------------------------------------
% by Bedrich Sousedik, October 2010.

bcdof = bc(1,:);
bcval = bc(2,:);

if length(bcval) ~= length(bcdof)
   error('Check assigned boundary conditions.')    
end 
% reorder (to have dofs in increasing/ascending order)
if ~issorted(bcdof) 
   [bcdof,index] = sort(bcdof,'ascend');
   bcval = bcval(index); 
end

ngdof = length(f);
bcdof_all = bcdof;

% take care of zero bc's first 
zdof = bcdof(bcval==0);
d = sparse(zdof,zdof,diag(A(zdof,zdof)),ngdof,ngdof);
A(zdof,:) = sparse(1,1,0,nnz(zdof),ngdof);
A(:,zdof) = sparse(1,1,0,ngdof,nnz(zdof));
A = A + d;
f(zdof) = zeros(nnz(zdof),1);
% DEBUG:
norm(sort(union(zdof,bcdof))-sort(bcdof_all));


% save for comparisons
%oldA = A;
%oldf = f;

% take care of non-zero b.c.'s next
bcdof = setdiff(bcdof,zdof);
bcval = bcval(ismember(bcdof_all,bcdof));
ndof = nnz(bcdof);
% extract diagonals
sc = diag(A(bcdof,bcdof));
% update RHS
f = f - A(:,bcdof)*bcval';
% remove rows and columns
A(bcdof,:) = spalloc(ndof,ngdof,1);
A(:,bcdof) = spalloc(ngdof,ndof,1);
% put values back on diagonal
A(bcdof,bcdof) = sparse(1:ndof,1:ndof,sc);
f(bcdof) = bcval'.*sc;

return

% old slow method
for i=1:length(bcdof)
   % get a scaling factor for the diagonal entries
   % (hopefully more numerically stable then simple 1 on the diagonal)
   sc = A(bcdof(i),bcdof(i));
   % 'remove' row
   A(bcdof(i),:) = spalloc(1,ngdof,1); %spalloc(ngdof,1,1);
   % 'remove' column and subtract it from the rhs
   %if bcval(i) ~= 0 - we took care of this above
      f = f - bcval(i)*A(:,bcdof(i));
   %end   
   A(:,bcdof(i)) = spalloc(ngdof,1,1);  %spalloc(1,ngdof,1);
   % put value on the diagonal 
   A(bcdof(i),bcdof(i)) = sc; %1;
   f(bcdof(i)) = bcval(i)*sc;
end
return % end of function

% % by Mon:
% n=length(bcdof);
% sdof=size(kk);
% for i=1:n
%     c=bcdof(i);
%     for j=1:sdof
%         ff(j)=ff(j)-bcval(i)*kk(j,c);
%         kk(c,j)=0;
%         kk(j,c)=0;
%     end
%     kk(c,c)=1;
%     ff(c)=bcval(i);
% end
