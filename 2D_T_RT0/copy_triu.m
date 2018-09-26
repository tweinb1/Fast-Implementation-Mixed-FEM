function A = copy_triu(A)

% COPY_TRIU
%   Copy the upper triangular part of given square matrices A such that
%   A(k,m,:) = A(m,k,:).
%
% SYNTAX:  A = copy_triu(A)
%
% IN:   A    matrices whose upper triangular parts we wish to copy
%
% OUT:  A    matrices with upper triangular parts copied to the lower
%            triangular parts (essentially making them symmetric)
%

nrows = size(A,1);
ncols = size(A,2);

if ( nrows ~= ncols )
    error('The given matrices are not square matrices.')
end

for m=1:nrows
    for k=m+1:nrows
        A(k,m,:) = A(m,k,:);
    end
end
