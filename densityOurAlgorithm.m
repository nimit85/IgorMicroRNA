% load data
load('4Haidar.mat', 'genes');

% z-score data = zero mean, unit variance
X = zscore(genes);
Q = X;

% parameters - neighborhood size and number of neighbors
epsil = 0.9;
minpts = 40;

%%

% adjacency is calculated based on epsilon distance and number of neighbors
adjList = zeros(size(Q,1),1);
reduceVal = 2;

for iter = 1 : size(Q,1)
    Y = Q;
    Y(iter,:) = [];
    cntr = 1;
    thisVal = sum(pdist2(Q(iter,:),Y) <= epsil);
    if thisVal >= minpts
        adjList(iter) = 1;
    else
        newpts = minpts;
        while ~adjList(iter)
            cntr = cntr + 1;
            newpts = floor(newpts/reduceVal);
            if ~newpts
                adjList(iter) = Inf;
                break;
            end
            if thisVal > newpts
                adjList(iter) = cntr;
            end
        end
    end
end

% determine what a constant "very large' is
S = sort(unique(adjList),'descend');
% S is number of clusters after algorithm
adjList(isinf(adjList)) = S(2)+1;
S(1) = S(2) + 1;

cntrOnS = zeros(size(S,1),1);
for j = 1 : size(S,1)
    cntrOnS(j) = sum(adjList == j);
end
minVal = min(cntrOnS);

% random sampling cntrOnS points from each cluster - samplesClust contains
% the randomly sampled points
samplesClust = zeros(size(S,1),minVal);
for j = 1 : size(S,1)
    F = find(adjList == j);
    onlyFew = F(randperm(size(F,1)));
    samplesClust(j,:) = onlyFew(1:minVal);
end

figure;
for j = 1 : max(adjList)
    plot(Q(adjList == j,1), Q(adjList == j,2), '.'); %[clrs(mod(j, length(clrs))+1) '.']);    
    hold all;
end
legend('show')
%%

C = B1H_0I_M0_R1_Jcel(samplesClust);
T = cell2table(C);
writetable(T,'ScatterPlotOutliers.txt','WriteRowNames',false,'WriteVariableNames',false);