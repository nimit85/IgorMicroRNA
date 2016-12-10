function intersectGenesAcrossMice
% get genes from the clustering algorithm and verify how many genes
% intersect across mice. 
% Task 1: Do the intersection across all 3 mice,
% Task 2: Do the intersection pair-wise (3 final results)
clear
clc
close all

% load cluster assignments
adjList1 = load(['adj_' num2str(1) '_List.mat']);
adjList1 = adjList1.adjList;
adjList2 = load(['adj_' num2str(2) '_List.mat']);
adjList2 = adjList2.adjList;
adjList3 = load(['adj_' num2str(3) '_List.mat']);
adjList3 = adjList3.adjList;

% load gene names
load('geneNames.mat');
geneNames = B1H_0I_M0_R1_Jcel;

% drop largest cluster, vertically concatenate rest of the genes
[UCluster1,idx1] = dropCluster(adjList1,geneNames);
[UCluster2,idx2] = dropCluster(adjList2,geneNames);
[UCluster3,idx3] = dropCluster(adjList3,geneNames);

% get the intersection of genes remaining after droppping largest cluster 
% from all mice
% interGenes = intersect(intersect(UCluster1,UCluster2),UCluster3);
% interIds = intersect(intersect(idx1,idx2),idx3);
interGenes = intersect(intersect(UCluster1,UCluster2),UCluster3);
writetable(table(interGenes),'finalGenesSubset.xlsx','WriteVariableNames',false);
interIds = intersect(intersect(idx1,idx2),idx3);

% assign the genes into signature classes based on expression levels,
% consider only the 4 conditions - Susceptible, Resistant, Susceptible +
% Resistant, Resistant + Susceptible
affyStruct1 = load('mouse1.mat');
F = fieldnames(affyStruct1); affyStruct1 = affyStruct1.(F{1});
affyStruct1 = affyStruct1(:,:,5:8);
affyStruct2 = load('mouse2.mat');
F = fieldnames(affyStruct2); affyStruct2 = affyStruct2.(F{1});
affyStruct2 = affyStruct2(:,:,5:8);
affyStruct3 = load('mouse3.mat');
F = fieldnames(affyStruct3); affyStruct3 = affyStruct3.(F{1});
affyStruct3 = affyStruct3(:,:,5:8);

% run through the gene set, if there is a majority vote for a signature
% capture that, if no majority vote, record those genes
% B, B with H, H, H with B
for jj = 1 : 4
    thisCondGenes1 = affyStruct1(interIds,:,jj);
    thisCondGenes2 = affyStruct2(interIds,:,jj);
    thisCondGenes3 = affyStruct3(interIds,:,jj);
    
    patt1 = takeDiffsCreateSign(thisCondGenes1);
    patt2 = takeDiffsCreateSign(thisCondGenes2);
    patt3 = takeDiffsCreateSign(thisCondGenes3);
    
    % do majority voting, or assign outlier label
    disp([patt1; patt2; patt3]);
    
end
end

function pattern = takeDiffsCreateSign(thisCondGene)
% find the difference in expression levels, use that to assign a class
% label
diff1 = thisCondGene(:,2) - thisCondGene(:,1);
diff2 = thisCondGene(:,3) - thisCondGene(:,2);
diff3 = thisCondGene(:,4) - thisCondGene(:,3);
pattern1 = sum(diff1) - 0.5*numel(diff1) > 0;
if pattern1 == 0
    pattern1 = -1;
end
pattern2 = sum(diff2) - 0.5*numel(diff2) > 0;
if pattern2 == 0
    pattern2 = -1;
end
pattern3 = sum(diff3) - 0.5*numel(diff3) > 0;
if pattern3 == 0
    pattern3 = -1;
end
pattern = [pattern1, pattern2, pattern3];
end

function [keepGenes, idxGenes] = dropCluster(adjList,geneNames)
finClusters = unique(adjList);
keepGenes = cell(0,0); idxGenes = [];
remCluster = 0;

% drop the largest cluster, for all other clusters - vertically concatenate
% the gene names and ids
for ii = 1 : max(finClusters)
    thisCluster = sum(adjList == ii);
    if remCluster < thisCluster        
        remCluster = thisCluster;
    else
        G = geneNames(adjList == ii);
        idx = find(adjList == ii);
        keepGenes = [keepGenes; G];
        idxGenes = [idxGenes; idx];
    end
end
end