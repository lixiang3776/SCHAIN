function [matrix,value,NcutDiscrete,t] = optimizeZ_PI(att_matrix1,mp_matrix,constraints,K,alpha,gamma,lambda,A)

P = size(mp_matrix,1);
att_num = size(att_matrix1,1);
n = size(mp_matrix,2);
matrix = zeros(n,n);

% lambda = ones(P,1);
% A = ones(att_num,1);

% alpha = 100000;
% beta = 10;
% gamma = 0.1; % gamma is regularization parameter

% for meta path similarity 
for i = 1:P
    matrix = matrix + lambda(i) * squeeze(mp_matrix(i,:,:));
end

matrix = alpha(1) * matrix;

for j = 1:att_num
    matrix = matrix + alpha(2) * A(j) * squeeze(att_matrix1(j,:,:));
end

matrix_new = matrix + constraints .* matrix;

% [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(matrix,K,constraints,beta);
% 
% clusters = zeros(n,1);
% 
% %get cluster index
% for i = 1:n
%     for j = 1:K
%         if NcutDiscrete(i,j) == 1
%             clusters(i) = j;
%             break;
%         end;
%     end;
% end;

% NcutDiscrete = clusters;

[NcutDiscrete,t] = PICk(matrix,matrix_new,K);
% [NcutDiscrete,t] = NJW(matrix,matrix_new,K);
% [NcutDiscrete,t] = NJW_IRAM(matrix,matrix_new,K);

clusters = NcutDiscrete;

sum_S = zeros(K,1);
sum_D = zeros(K,1);
D = sum(matrix,2);

for i = 1:n
    sum_D(clusters(i)) = sum_D(clusters(i)) + D(i);
    for j = 1:n
        if clusters(i) == clusters(j)
           sum_S(clusters(i)) = sum_S(clusters(i)) + matrix(i,j);
           if constraints(i,j) == 1
               sum_S(clusters(i)) = sum_S(clusters(i)) + matrix(i,j);
           elseif constraints(i,j) == -1
               sum_S(clusters(i)) = sum_S(clusters(i)) - matrix(i,j);
           end;
        end;
    end;
end;

value = K - (sum(sum_S./sum_D) - gamma * (norm(lambda)^2 + norm(A)^2));

disp(sum(NcutDiscrete));

disp(value);

% [lambda, A] =  optimizeWeights(att_matrix,mp_matrix,NcutDiscrete,K,constraints,alpha,beta,gamma);

end