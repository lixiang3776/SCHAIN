function [nmi] = calculateNMI_semi(cluster,true_cluster,seeds)
n = size(cluster,1);
K = size(cluster,2);
I = 0;
for i=1:K    %cluster
    for j=1:K    %label
        wkcj = 0;
		wk = 0;
		cj = 0;
		N = 0;
        for p=1:n
            if ~ismember(p,seeds)
                N = N + 1;
                if (cluster(p,i) == 1) && (true_cluster(p) == j)
                   wkcj = wkcj + 1;
                end;
                if cluster(p,i) == 1
                   wk = wk + 1;
                end;
                if true_cluster(p) == j
                   cj = cj + 1;
                end;
            end;
        end;
        if wkcj > 0
			I = I + (wkcj/N)*log((N*wkcj)/(wk*cj) + realmin);
        end;
    end;
end;

H1 = 0;
for i=1:K
    wk = 0;
    N = 0;
    for p=1:n
        if ~ismember(p,seeds)
            N = N + 1;
            if cluster(p,i) == 1
                wk = wk + 1;
            end;
        end;
    end;
    if wk > 0
        H1 = H1 + (wk/N)*log(wk/N + realmin);
    end;
end;
H1 = H1 * (-1);

H2 = 0;
for j=1:K
    cj = 0;
    N = 0;
    for p=1:n
        if ~ismember(p,seeds)
            N = N + 1;
            if true_cluster(p) == j
                cj = cj + 1;
            end;
        end;
    end;
    if cj > 0
        H2 = H2 + (cj/N) * log(cj/N + realmin);
    end;
end;
H2 = H2 * (-1);
nmi = 2 * I / (H1 + H2);
% fprintf('%f %f %f', I, H1, H2);
end

