% A version of the Euclidean k-means algorithm where the several trials are
% ran and the most frequent result is chosen as the final clustering
%
% Input:
% data - row-instance data matrix
% k - number of clusters
% trials - number of trials to run before determining best clustering
% implementation - which specific k-means implementation to use
%
% implementation:
% 'm' - one that comes with matlab
% 'cl' - one by Chen & Lin (wychen@alumni.cs.ucsb.edu)
%
% Output:
% pred - predicted cluster labels
%
% Author: Frank Lin (frank@cs.cmu.edu)

function pred=kmeans_freq(data,k,trials,implementation)

fprintf('running %d k-means trials\n',trials);

% data size
n=size(data,1);

% initialize prediction store
preds=zeros(n,trials);

for i=1:trials
    
    % run k-means
    if strcmp(implementation,'m')
        % MATLAB version
        IDX=kmeans(data,k,'emptyaction','singleton');
    elseif strcmp(implementation,'cl')
        % Chen & Lin version
        IDX=kmeans_cl(data,'random',k);
    else
        fprintf('implementation not recognized: %s\n',implementation)
    end

    % normalize labels
    preds(:,i)=normlabels(IDX);
    
end

% find unique clusterings and their indices
[upreds,~,uinds]=unique(preds','rows');

% find clustering counts
predc=histc(uinds,1:size(upreds,1));

% find most frequent index and count
[freqc,freqi]=max(predc);

fprintf('most frequent clustering: %d/%d (%f)\n',freqc,trials,freqc/trials);

pred=upreds(freqi,:)';

end
