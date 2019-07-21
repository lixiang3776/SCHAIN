function [lambda, A, mu_value] = optimizeWeights_testPI(att_matrix1, mp_matrix, clusters, K, constraints, alpha, gamma)
% % The function is used to learn weights of meta paths and attributes.
% Input:
% att_matrix: n*att_num, where att_num is the number of attributes.
% mp_matrix: P*n*n, where P is the number of meta paths used and n is the
% number of objects. Each mp_matrix(i,:,:) is a meta path based
% similarity matrix.
% NcutDiscrete: n*1, is the clustering result of z in the previous step. 
% K: the number of clusters.
% constraints: the constraint matrix of n*n.
% alpha: the balance parameter.
% gamma: the regularization parameter.
% 
% Output:
% lambda: the meta path parameter vector.
% A: the attribute paramter vector.
% mu_value: the value of mu.
% 
% % 

P = size(mp_matrix,1);
att_num = size(att_matrix1,1);
para_num = P+att_num;
n = size(mp_matrix,2);

X = sym('x',[1 para_num]);
assume(X, 'real');

mp_sum = zeros(K,P);% K clusters, P meta paths
mp_sum_D = zeros(K,P);% K clusters, P meta paths
attr_sum = zeros(K,att_num);% K clusters, att_num attributes
attr_sum_D = zeros(K,att_num);% K clusters, att_num attributes


% for meta paths
for p = 1:P
    mp = squeeze(mp_matrix(p,:,:));
    mp = alpha(1) * mp;
    for i = 1:n
        for j = i+1:n
            if clusters(i) == clusters(j)
                mp_sum(clusters(i),p) = mp_sum(clusters(i),p) + mp(i,j);
                mp_sum(clusters(j),p) = mp_sum(clusters(j),p) + mp(i,j);
                if constraints(i,j) == 1
                    mp_sum(clusters(i),p) = mp_sum(clusters(i),p) + mp(i,j);
                    mp_sum(clusters(j),p) = mp_sum(clusters(j),p) + mp(i,j);
                elseif constraints(i,j) == -1
                    mp_sum(clusters(i),p) = mp_sum(clusters(i),p) - mp(i,j);
                    mp_sum(clusters(j),p) = mp_sum(clusters(j),p) - mp(i,j);
                end;
            end;
            mp_sum_D(clusters(i),p) = mp_sum_D(clusters(i),p) + mp(i,j);
            mp_sum_D(clusters(j),p) = mp_sum_D(clusters(j),p) + mp(i,j);
        end;
    end;
end;

% for attributes
for q = 1:att_num
    att = squeeze(att_matrix1(q,:,:));
    att = alpha(2) * att;
    for i = 1:n
        for j = i+1:n
            if clusters(i) == clusters(j)
                attr_sum(clusters(i),q) = attr_sum(clusters(i),q) + att(i,j);
                attr_sum(clusters(j),q) = attr_sum(clusters(j),q) + att(i,j);
                if constraints(i,j) == 1
                    attr_sum(clusters(i),q) = attr_sum(clusters(i),q) + att(i,j);
                    attr_sum(clusters(j),q) = attr_sum(clusters(j),q) + att(i,j);
                elseif constraints(i,j) == -1
                    attr_sum(clusters(i),q) = attr_sum(clusters(i),q) - att(i,j);
                    attr_sum(clusters(j),q) = attr_sum(clusters(j),q) - att(i,j);
                end;
            end;
            attr_sum_D(clusters(i),q) = attr_sum_D(clusters(i),q) + att(i,j);
            attr_sum_D(clusters(j),q) = attr_sum_D(clusters(j),q) + att(i,j);
        end;
    end;
end;

Q = 0;

for i=1:K
    num = 0;
    den = 0;
    for p = 1 : P
        num = num + sym(sprintf('%d',mp_sum(i,p))) * X(p);
    end;
    for j = 1 : att_num
        num = num + sym(sprintf('%d',attr_sum(i,j))) * X(P+j);
    end;
    for p = 1 : P
        den = den + sym(sprintf('%d',mp_sum_D(i,p))) * X(p);
    end;
    for j = 1 : att_num
        den = den + sym(sprintf('%d',attr_sum_D(i,j))) * X(P+j);
    end;

    Q = Q + num/den;
end;

for h = 1 : para_num
    Q = Q - sym(sprintf('%f',gamma)) * X(h)^2;
end;

Q = Q + 2 * gamma;

[num(X), den(X)] = numden(Q);

numerator(X) = expand(num);
denominator(X) = expand(den);

global g0;
global g1;
g0 = matlabFunction(numerator);
g1 = matlabFunction(denominator);

global mu;
mu  = 0;
pre_mu = realmax; % pre_mu records the previous mu.
max_mu = realmax; % max_mu is the maximum mu leading to a negative value of the function.
lambda = ones(P,1)./P;
A = ones(att_num,1)./att_num;


while 1
    if abs(mu-pre_mu) < 1e-5
        break;
    end;
    
    if mu < max_mu
        x0 = rand(1,para_num);
        temp1 = [ones(1,P),zeros(1,att_num)];
        temp2 = [zeros(1,P),ones(1,att_num)];
        Aeq = [temp1;temp2];
        beq = [1;1];
        lb = zeros(1,para_num);
        ub = ones(1,para_num);
        nonlcon = [];
        options = optimoptions('fmincon','MaxFunctionEvaluations',8000,'ObjectiveLimit',-1e40,'Display','off');
        options.Algorithm = 'sqp';
        [x, fval, exitflag] = fmincon(@fun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
        temp_value = fval*(-1); % temp_value records the value of the function. 
        
        if (exitflag == 1) && (temp_value >= 0)% if the value >= 0, it indicates that we need to further increase mu
            pre_mu = mu;
            lambda = x(1:P)';
            A = x(P+1:end)';     
            arg_x = {};
            for len = 1:length(x)
                arg_x = [arg_x,x(len)];
            end
            mu = double(numerator(arg_x{:})/denominator(arg_x{:}));
            disp(exitflag);
            disp(pre_mu);
            disp(mu);
        elseif (exitflag == 1) && (temp_value < 0)% if the value < 0, it indicates that we need to further decrease mu
            disp(exitflag);
            disp(temp_value);        
            if max_mu > mu
                max_mu = mu;
            end;
            mu = (pre_mu + max_mu)/2;
            disp(mu);
            disp(max_mu);
        else % in other cases, just neglect the result
%             disp(x);
%             disp(exitflag);
%             disp(temp_value);
%             disp(mu);
        end;
    else % if mu >= max_mu, we directly discard it and derive a new smaller mu
        mu = (pre_mu + max_mu)/2;
    end;
end;

mu_value = mu;

end

function f = fun(y)
x = {};
for i=1:length(y)
    x = [x, y(i)];
end;

global g0;
global g1;
global mu;

f = mu * g1(x{:}) - g0(x{:});
end

