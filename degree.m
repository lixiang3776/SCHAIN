% Creates a diagonal degree matrix D from W, where D(i,i) is the 
% row-sum of the i-th row of W. D will be a sparse matrix.
%
% Author: Frank Lin (frank@cs.cmu.edu)

function D=degree(W)

n=size(W,1);
D=spdiags(sum(W,2),0,n,n);

end