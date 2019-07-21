% Power Iteration Clustering projection
%
% Input:
% W - row-normalized affinity matrix
% v0 - starting vector
% conv - convergence threshold
% maxit - maximum number of iterations
%
% Output:
% vt - 1-d PIC embedding
% i - iterations ran
% t - runtime
%
% Author: Frank Lin (frank@cs.cmu.edu)

function [vt,i,t]=pic(W,v0,conv,maxit)

n=size(W,1);

vt=v0;
dt=ones(n,1);
dtp=zeros(n,1);

i=0;

% tic;
while(max(abs(dt-dtp))>conv&&i<maxit)
    
    vtp=vt;
    dtp=dt;
    
    vt=W*vt;
    vt=vt/sum(vt);
    
    dt=abs(vt-vtp);
    
    i=i+1;
    
end
% t=toc;

end