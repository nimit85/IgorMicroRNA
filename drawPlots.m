folderNames = {'MTB_B','MTB_H','MTB+IFN_B','MTB+IFN_H'};
mice = {'R1','R2','R3'};

clusterFolder = 'output/';
plotsFolder = 'plots/';
clustMice1 = dir([clusterFolder 'mouse_1_cluster_*.xlsx']);
clustMice2 = dir([clusterFolder 'mouse_2_cluster_*.xlsx']);
clustMice3 = dir([clusterFolder 'mouse_3_cluster_*.xlsx']);

xAxisLabels = [0,6,30,78];

for ii = 1 : size(folderNames,2)
   for jj = 1 : size(mice,2)
       % load text file with probe information
       fileName = [folderNames{ii} '_' mice{jj} '.txt'];
       probes = readtable([folderNames{ii},'/',mice{jj},'/',fileName]);
       
       if jj == 1
           miceFolder = clustMice1;
       elseif jj == 2
           miceFolder = clustMice2;
       else
           miceFolder = clustMice3;
       end
       
       % load cluster files
       for kk = 1 : size(miceFolder,1)
           clusts = readtable([clusterFolder,miceFolder(kk).name]);
           
           % find all bin assignments
           U = unique(clusts.bins);
           
           for ll = 1 : size(U,1)
           
               % find intersection between cluster assignments and probes to
               % find the time point information
               X1 = probes(:,1);
               thisIdx = clusts.bins == U(ll);
               X2 = clusts(thisIdx,1);
               
               [~,ia,~] = intersect(X1.Var1,X2.Name);
               timePts = table2array(probes(ia,2:5));

               % get probe names
               legEntries = table2array(probes(ia,1));

               fig = figure('Visible','off');
               plot(xAxisLabels,timePts);
               ax = gca;
               ax.XTick = xAxisLabels;
               legend(legEntries,'Location','EastOutside');
               newName = strrep(miceFolder(kk).name,'.xlsx','');
               newName = [newName '_bins_' num2str(U(ll)) '.jpeg'];
               print(fig,'-djpeg',[plotsFolder,mice{jj},'/',folderNames{ii},'/',newName]);
           end
       end
   end
end