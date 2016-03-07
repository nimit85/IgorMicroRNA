% load data
load('4Haidar.mat', 'genes');

% z-score data = zero mean, unit variance
X = zscore(genes);
Q = X;

% parameters - neighborhood size and number of neighbors
epsil = 0.1;
minpts = 1;

%%

% set parameters for DBSCAN, run DBSCAN

% we dont need this part anymore
[idx, scores] = dbscan(Q, epsil, minpts);

%%

% plot results of DBSCAN
figure; hold on;
for j = 1:max(idx)
    plot(Q(idx==j,1), Q(idx==j,2), '.'); %[clrs(mod(j, length(clrs))+1) '.']);    
end

% determine what a constant "very large' is
S = sort(unique(idx),'descend');

cntrOnS = zeros(size(S,1),1);
for j = 1 : size(S,1)
    cntrOnS(j) = sum(idx == j);
end
minVal = min(cntrOnS);

% random sampling cntrOnS points from each cluster - samplesClust contains
% the randomly sampled points
samplesClust = zeros(size(S,1),minVal);
for j = 1 : size(S,1)
    F = find(idx == j);
    onlyFew = F(randperm(size(F,1)));
    samplesClust(j,:) = onlyFew(1:minVal);
end