function [idx,t,E] = PICk(matrix,matrix_new,k)
t_start = tic;
baseconv=1e-5;
% n=size(matrix_new,1);
% W=normrow(matrix_new);

Dinv=degree(matrix);
n = size(Dinv,1);
for i=1:n
  Dinv(i,i)=1/sqrt(Dinv(i,i));
end
W=Dinv*matrix_new*Dinv;
E=zeros(k,n);

for i=1:k
    v0=rand(n,1);
    [e,~]=pic(W,v0,baseconv/n,1000);
    E(i,:)=e;
end
t = toc(t_start);
[idx] = postprocess(E',k);

end

