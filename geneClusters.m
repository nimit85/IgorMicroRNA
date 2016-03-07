clc; close all; clear

% Read gene file
fileId = fopen('M1_Gene_Loadings_Data.txt','r');
genes = textscan(fileId,'%u %f %f');
genes = [genes{:,2}, genes{:,3}];
%%

% Maintain random seeds for reproducibility
rng default;

% Clustering algorithm
K = 10:40;
opts = statset('Display','final','MaxIter',10000);

indices = zeros(size(genes,1),length(K));
for iter = K
    indices(:,iter-K(1)+1) = kmeans(genes,iter,'Distance','cityblock','Options',opts,...
        'Replicates',10);
end

%%

% Evaluate clustering solutions to identify optimal K value
eva = evalclusters(genes, indices, 'CalinskiHarabasz');

%% 

% Silhouettes plots to estimate best number of clusters
silh = zeros(size(genes,1),length(K));
for iter = K
    fig = figure('Visible','off');
    [silh(:,iter-K(1)+1), ~] = silhouette(genes,indices(:,iter-K(1)+1),'cityblock');
    h = gca;
    h.Children.EdgeColor = [.8 .8 1];
    xlabel('Silhouette Value'); ylabel('Cluster');
    print(fig,'-r300','-dpng',['cluster_with_' num2str(iter) '.png']);
    disp(iter);
end

%%

% Scatter plots
for iter = K
    fig = figure('visible','off');
    gscatter(genes(:,1), genes(:,2),indices(:,iter-3));
    print(fig,'-r300','-dpng',['scatter_with_' num2str(iter) '.png']);
end