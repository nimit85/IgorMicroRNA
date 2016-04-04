% load outliers file
someGenes = load('HH_BB1_gene_outliers.txt');
load('geneNames.mat');

% choose some number of genes
K = 20;

% sort first and second columns
[col1, idx1] = sort(someGenes(:,1));
[col2, idx2] = sort(someGenes(:,2));

% get gene names
C2 = B1H_0I_M0_R1_Jcel(idx2(1:K));
C1 = B1H_0I_M0_R1_Jcel(idx1(1:K));

% sum of the 2 columns
col3 = someGenes(:,1) + someGenes(:,2);

% sort the sum, get gene names
[col3, idx3] = sort(col3);
C3 = B1H_0I_M0_R1_Jcel(idx3(1:K));
