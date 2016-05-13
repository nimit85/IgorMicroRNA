function corrAnalysis(affyStruct1,opts)
close all

% load cluster assignments
load adjList.mat

affyStruct1 = reshape(affyStruct1,size(affyStruct1,1),...
    size(affyStruct1,2)*size(affyStruct1,3));

% load the dimensionality reduced data
load('HH_BB1_gene_loadings.mat', 'genes');

% load gene names
load('geneNames.mat');
geneNames = B1H_0I_M0_R1_Jcel;

% drop largest cluster, find correlation coefficients for rest of clusters
U = unique(adjList);

remCluster = 0;
markDel = 0;
for ii = 1 : max(U)
    thisCluster = sum(adjList == ii);
    if remCluster < thisCluster
        markDel = ii;
        remCluster = thisCluster;
    end
end

for ii = 1 : max(U)
    if ii ~= markDel
        % get adjacency list for this cluster
        thisCluster = adjList == ii;
        corrStruct = affyStruct1(thisCluster,:);
        
        % calculate correlation coefficient for this cluster
        corrMat = triu(corr(corrStruct','type','Pearson'));
        
        % make diagonal entries 0
        corrMat(logical(eye(size(corrMat)))) = 0;
        
        % get signed version of correlation matrix
        simMat = corrMat; %0.5*(1 + corrMat);
        
        % hard threshold on correlation
        params = 0.65:0.05:0.95;
        % soft threshold on adjacency function
        betaVal = 1:20;
        slope = zeros(length(params),length(betaVal));
        meanConn = zeros(length(params),length(betaVal));
        Rsquared = zeros(length(params),length(betaVal));
        
        for jj = 1 : length(params)
            for kk = 1 : length(betaVal)
                adjMat = simMat .^ betaVal(kk);
                [Rrow, Rcol] = find(adjMat >= params(jj));
                
                % get the degree of nodes in the graph
                src = Rrow;
                tgt = Rcol;
                graphG = graph(src,tgt);
                degs = degree(graphG);
                
                % see if the degree distribution follows a power law
                U = unique(degs);
                if isempty(U)
                    continue;
                end
                [Hcounts, ~] = histcounts(degs,U);
                Hcounts = Hcounts/sum(Hcounts);
                
                y = log10(U(2:end)');
                x = log10(Hcounts);
                % slope and intercept for the log-log fit
                regressParams = polyfit(x,y,1);
                
                % calculate R^2 statistic, mean connectivity
                yHat = polyval(regressParams, x);
                SStot = norm(y - mean(y))^2;
                SSres = norm(y - yHat)^2;
                Rsquared(jj,kk) = 1 - SSres/SStot;
                meanConn(jj,kk) = mean(degs);
                slope(jj,kk) = regressParams(1);
            end
            
            
            % [beta(jj), ~, loglike(jj)] = plfit(Hcounts,'range',1.5:0.05:3);
        end
        
        % choose best combination - R^2 > 0.8, high mean connectivity,
        % slope about -1
        keepR = Rsquared > 0.8;
        keepHigh = meanConn > 1;
        keepSlope = abs(slope) > 0.7 & abs(slope) < 1.3;
        bestCombo = and(and(keepR,keepHigh),keepSlope);
        
        [rownum, colnum] = find(bestCombo == 1);
        bestParams = params(rownum);
        bestBeta = betaVal(colnum);
        
        for jj = 1 : length(bestParams)
            adjMat = simMat .^ bestBeta(jj);
            [Rrow, Rcol] = find(adjMat >= bestParams(jj));
            
            % create gene graph
            src = Rrow;
            tgt = Rcol;
            
            graphG = graph(src,tgt);
            genesCluster = geneNames(thisCluster);
            F = degree(graphG);
            
            % remove nodes with degree 0 aka isolated nodes
            remNodes = find(F == 0);
            genesCluster(remNodes) = [];
            graphG = rmnode(graphG,remNodes);
            
            % get gene names for nodes
            srcNodes = graphG.Edges.EndNodes(:,1);
            tgtNodes = graphG.Edges.EndNodes(:,2);
            srcNames = genesCluster(srcNodes);
            tgtNames = genesCluster(tgtNodes);
            unqNames = unique([srcNames; tgtNames]);
            graphG.Nodes = unqNames;
            
            % to get paths in graph
            bins = conncomp(graphG,'OutputForm','cell');
            save(['mouse_' num2str(opts) '_genebins_' num2str(U) ...
                '_params_' num2str(jj)],'bins');
            
            % plot gene graph and heat map
            fig1 = figure('Visible','off'); 
            plot(graphG,'NodeLabel',unqNames,'NodeLabelMode','auto');
            title(['Params: ' num2str(bestBeta(jj)) ' ' num2str(bestParams(jj))]);
            savefig(fig1,['mouse_' num2str(opts) '_cluster_' num2str(U) ...
                '_params_' num2str(jj)]);
            fig2 = figure('Visible','off'); 
            imagesc(adjMat); axis equal; axis tight; colormap(jet);
            title(['Params: ' num2str(bestBeta(jj)) ' ' num2str(bestParams(jj))]);
            savefig(fig2,['mouse_' num2str(opts) '_heatmap_cluster_' num2str(U) ...
                '_params_' num2str(jj)]);
        end        
    end
end
end