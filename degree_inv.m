% Creates a diagonal degree matrix D from W, where D(i,i) is the reciprocal
% of the row-sum of the i-th row of W. D will be a sparse matrix.
%
% Author: Frank Lin (frank@cs.cmu.edu)

function D=degree_inv(W)

n=size(W,1);
d=sum(W,2).^-1;
d(isinf(d))=0;
D=spdiags(d,0,n,n);

end