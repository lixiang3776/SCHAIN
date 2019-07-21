function [final_nmi,origin_nmi,NcutDiscrete] =  optimization_testPI(att_matrix1, mp_matrix, K,true_cluster,percentage)
tic;
out_it = 0;
gamma = 0.5;

P = size(mp_matrix,1);
att_num = size(att_matrix1,1);
n = size(att_matrix1,2);

lambda = ones(P,1)./P;
A = ones(att_num,1)./att_num;
final_nmi = 0;
origin_nmi = 0;
% final_accuracy = 0;
samples_num = 10;
time_ncuts = 0;
num_ncuts = 0;
time_weights = 0;
num_weights = 0;

while out_it < samples_num
%     it = 0;
%     lambda = ones(P,1)./P;
%     A = ones(att_num,1)./att_num;
    constraints = zeros(n,n); 
    num = floor(n * percentage);
    seeds = randperm(n,num);
    for i=1:num
        for j=i+1:num
            if true_cluster(seeds(i)) == true_cluster(seeds(j))
                constraints(seeds(i),seeds(j)) = 1;
                constraints(seeds(j),seeds(i)) = 1;
            else
                constraints(seeds(i),seeds(j)) = -1;
                constraints(seeds(j),seeds(i)) = -1;
            end;
        end;
    end;

    nmi = 0;
    value_max = realmax;
    it = 0;
    lambda = ones(P,1)./P;
    A = ones(att_num,1)./att_num;
    alpha = [0.5;0.5];
    ITERATION = 5;
    while it < ITERATION
        %     && epsi > EPSI_T
        tStart1 = tic;
        [matrix,value,NcutDiscrete,t] =  optimizeZ_PI(att_matrix1,mp_matrix,constraints,K,alpha,gamma,lambda,A);
        
        tElapsed1 = toc(tStart1);
        
        if value >= value_max
            break;
        else
            value_max = value;
        end
        
        time_ncuts = time_ncuts + t;
        num_ncuts = num_ncuts + 1;
              
        nmi_temp = calculateNMI_semiPI(NcutDiscrete, true_cluster,seeds)
        %     accuracy_temp = calculateAccuracy(NcutDiscrete,true_cluster)
        if it == 0
            origin_nmi = origin_nmi + nmi_temp;
        end;
        
        if it==0 && out_it < 1
            fileID = fopen('exp.txt','w');
            fprintf(fileID,'%d nmi: %f\n',(it+1),nmi_temp);
            %         fprintf(fileID,'%d accuracy: %f\n',(it+1),accuracy_temp);
            fprintf(fileID,'%d value: %f\n',(it+1),value);
            fprintf(fileID,'\n');
            fclose(fileID);
        else
            fileID = fopen('exp.txt','a');
            fprintf(fileID,'%d nmi: %f\n',(it+1),nmi_temp);
            %         fprintf(fileID,'%d accuracy: %f\n',(it+1),accuracy_temp);
            fprintf(fileID,'%d value: %f\n',(it+1),value);
            fprintf(fileID,'\n');
            fclose(fileID);
        end;
        
        if it < (ITERATION-1)
            
            tStart2 = tic;
            [lambda,A] =  optimizeWeights_testPI(att_matrix1, mp_matrix, NcutDiscrete, K, constraints, alpha, gamma);
            tElapsed2 = toc(tStart2);
            
            time_weights = time_weights + tElapsed2;
            num_weights = num_weights + 1;
            
%             fprintf('optimize alpha:');
%             [alpha, mu_value] = optimizeWeights_alphaPI(att_matrix1, mp_matrix, NcutDiscrete, K, constraints, lambda, A, gamma);
            
            fileID = fopen('exp.txt','a');
            fprintf(fileID,'lambda: ');
            fprintf(fileID,'%f ',lambda);
            fprintf(fileID,'\n');
            fprintf(fileID,'A: ');
            fprintf(fileID,'%f ',A);
            fprintf(fileID,'\n');
%             fprintf(fileID,'alpha: ');
%             fprintf(fileID,' %f',alpha);
%             fprintf(fileID,'%f %f',fun_value,value1);
%             fprintf(fileID,' %f %f',pre_mu,pre_value);
%             fprintf(fileID,' %f %f',mu,max_mu);
            fprintf(fileID,'\n');
            fclose(fileID);
        end;
        
%         if it == (ITERATION - 1)
            nmi = nmi_temp;
            %         accuracy = accuracy_temp;
%         end;
        it = it + 1;
    end;
    
    final_nmi = final_nmi + nmi;
    % final_accuracy = final_accuracy + accuracy;
    out_it = out_it + 1;
end;
toc;
final_nmi = final_nmi/samples_num;
origin_nmi = origin_nmi/samples_num;

fprintf('time for ncuts: %f, number: %d\n', time_ncuts/num_ncuts, num_ncuts);
fprintf('time for weight learning: %f, number: %d\n', time_weights/num_weights, num_weights);

fileID = fopen('exp.txt','a');
fprintf(fileID,'%f : %f\n',final_nmi,origin_nmi);
fclose(fileID);
end