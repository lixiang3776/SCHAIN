function [idx,t_end] = NJW(matrix,matrix_new,k)

% W=normlapsym(full(A));
% [E,V]=eig(W);
% [~,index]=sort(diag(V));
% E=E(:,index);
% E=E(:,1:k);
t_start = tic;
Dinv=degree(matrix);
n = size(Dinv,1);
for i=1:n
  Dinv(i,i)=1/sqrt(Dinv(i,i));
end

N=Dinv*matrix_new*Dinv;
[E,V]=eig(N); % eigendecomposition on N and choose the largest eigenvectors
[~,index]=sort(-diag(V));
E=E(:,index);
E=E(:,1:k);
E = Dinv * E;
t_end=toc(t_start);
[idx] = postprocess(E,k);



