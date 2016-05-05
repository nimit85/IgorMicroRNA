load('mouse3.mat')
%BigY=dataset(affyStruct3);
Y=dataset(affyStruct3);
clear affyStruct3
 %Y= BigY(:,:,5:8);

label2 = {'0'
    '6'
    '34'
    '78'
    };
Y.label{2} = label2;

label3 = {
 %'MTB+IFN_B'
%'MTB+IFN_B1H'
%'MTB+IFN_H'
%'MTB+IFN_H1B'
'MTB_B'
'MTB_B1H'
'MTB_H'
'MTB_H1B'
    };



Y.label{3} = label3;
Y.title{1}= 'Genes'
Y.title{3}= 'Experiments';
Y.title{2}= 'Time';

editds(Y);
Fac = 2; 
[Factors] = parafac(Y,Fac);  %


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get gene loadings
genes= Factors.loads{1};


%[x indices] = sort(A);



% Maintain random seeds for reproducibility
rng default;

% Clustering algorithm
K = 3:8;

opts = statset('Display','final','MaxIter',10000);

indices = zeros(size(genes,1),length(K));
for iter = K
    indices(:,iter-K(1)+1) = kmeans(genes,iter,'Distance','cityblock','Options',opts,...
        'Replicates',10);
end
%

% Silhouettes plots to estimate best number of clusters
silh = zeros(size(genes,1),length(K));
for iter = K
    %fig = figure('Visible','off');
    [silh(:,iter-K(1)+1), ~] = silhouette(genes,indices(:,iter-K(1)+1),'cityblock');
    %h = gca;
    %h.Children.EdgeColor = [.8 .8 1];
    %xlabel('Silhouette Value'); ylabel('Cluster');
    %print(fig,'-r300','-dpng',['cluster_with_' num2str(iter) '.png']);
    disp(iter);
end
means = zeros(1:6);
for i = 6
k_means(i) = mean(silh(:,i));
end

%%

% Scatter plots
for iter = K
    fig = figure('visible','off');
    gscatter(genes(:,1), genes(:,2),indices(:,iter-K(1)+1));
    print(fig,'-r300','-dpng',['scatter_with_' num2str(iter) '.png']);
end
means = zeros(1:21);
for i = K
means(i) = mean(silh(:,i));
end

