function [idx,t2] = NJW_IRAM(matrix,matrix_new,k)

% W=normlapsym(full(matrix));
% [E,V]=eig(W);
% [~,index]=sort(diag(V));
% E=E(:,index);
% E=E(:,1:k);

% W=normlapsym(matrix);
% opts.disp=0;
% [E,~]=eigs(W,k,'SM',opts);
t1 = tic;
Dinv=degree(matrix);
n = size(Dinv,1);
for i=1:n
  Dinv(i,i)=1/sqrt(Dinv(i,i));
end

N=Dinv*matrix_new*Dinv;
[E,~]=eigs(N,k,'LM');
E = Dinv * E;
t2 = toc(t1);
[idx] = postprocess(E,k);




