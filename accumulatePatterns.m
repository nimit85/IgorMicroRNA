folderNames = {'MTB_B','MTB_H','MTB+IFN_B','MTB+IFN_H'};
mice = {'R1','R2','R3'};

clusterFolder = 'output/';
resultsFolder = 'results/';
clustMice1 = dir([clusterFolder 'mouse_1_cluster_*.xlsx']);
clustMice2 = dir([clusterFolder 'mouse_2_cluster_*.xlsx']);
clustMice3 = dir([clusterFolder 'mouse_3_cluster_*.xlsx']);

% for all mice
for jj = 1 : size(mice,2)
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
            % plot the accumulated pattern
            fig = figure('Visible','off'); hold on;
    
            % find intersection between cluster assignments and probes to
            % find the time point information
            % for all 4 conditions
            cntr = 0; toSaveFlag = 1;
            for ii = 1 : size(folderNames,2)
                
                % load text file with probe information
                fileName = [folderNames{ii} '_' mice{jj} '.txt'];
                probes = readtable([folderNames{ii},'/',mice{jj},'/',fileName]);
                
                X1 = probes(:,1);
                thisIdx = clusts.bins == U(ll);
                X2 = clusts(thisIdx,1);
                
                [~,ia,~] = intersect(X1.Var1,X2.Name);
                timePts = table2array(probes(ia,2:5));
                
                if size(timePts,1) < 10
                    toSaveFlag = 0;
                    break;
                end
                
                % create accumulated pattern
                diff1 = (timePts(:,2) - timePts(:,1)) > 0;
                diff2 = (timePts(:,3) - timePts(:,2)) > 0;
                diff3 = (timePts(:,4) - timePts(:,3)) > 0;
                
                pattern1 = sum(diff1) > 0.5*numel(diff1);
                pattern2 = sum(diff2) > 0.5*numel(diff2);
                pattern3 = sum(diff3) > 0.5*numel(diff3);
                
                pattern = [pattern1, pattern2, pattern3];
                timePts = [1,2,3];
                
                % get probe names
                legEntries = probes(ia,1);
                plot(timePts,pattern+cntr);
                cntr = cntr + 2;
            end
            
            if toSaveFlag
                xlabel('Time difference patterns');
                ylabel('Accumulated pattern');
                ax = gca;
                ax.XTickLabel = [1 2 3];
                legend(folderNames,'Location','EastOutside');
                ax.YTickLabel = [0 1 0 1 0 1 0 1];
                newName = strrep(miceFolder(kk).name,'.xlsx','');
                newName = [newName '_bins_' num2str(U(ll)) '.jpeg'];
                print(fig,'-djpeg',[resultsFolder,mice{jj},'/',newName]);
                writetable(legEntries,[resultsFolder,mice{jj},'/',...
                    strrep(newName,'.jpeg','.xlsx')],...
                    'WriteVariableNames',false,'WriteRowNames',false);                
            end
            hold off; close all;
        end
    end
end