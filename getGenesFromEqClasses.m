srcDir = 'results/'; % source directory
fid = fopen([srcDir, 'patternNames.txt']); % load the equivalence class file
C = textscan(fid,'%s %s %s -- %s'); % throw away first read

genes = [];
eqClass = 1;

while ~feof(fid) % exit when full file is read
    disp(eqClass);
    C = textscan(fid,'%s %s %s -- %s');
    if strcmpi(C{1},'Equivalence')
        if ~isempty(genes)
            % check to see which genes appear in majority of clusters
            
        end
        continue;
    else
        % load all the files, vertically concatenate genes
        fNames = C{1};
        fNames = fNames(1:end-1);
        for ii = 1 : size(fNames,1)
            
            % check which mouse folder the file is coming from and load the
            % file accordingly
            if strfind(fNames{ii},'mouse_1')
                X = readtable([srcDir, 'R1/', fNames{ii}], ...
                    'ReadVariableNames',false);
                genes = [genes; X];
            elseif strfind(fNames{ii},'mouse_2')
                X = readtable([srcDir, 'R2/', fNames{ii}], ...
                    'ReadVariableNames',false);
                genes = [genes; X];
            else
                X = readtable([srcDir, 'R3/', fNames{ii}], ...
                    'ReadVariableNames',false);
                genes = [genes; X];
            end
        end
        
        % convert to categorical array
        genes = categorical(genes.Var1);
        
        % find unique values
        U = unique(genes);
        
        % run histogram-type code to see which genes occur frequently
        keepGenes = [];
        for ii = 1 : size(U,1)
            idx = find(U(ii) == genes);
            if size(fNames,1) < 10
                % for small equivalence classes, a larger threshold works
                if numel(idx) >= 0.5*size(fNames,1) % 50% of files
                    keepGenes = [keepGenes; U(ii)];
                end
            else
                % for larger equivalence classes, a smaller threshold is
                % required.
                if numel(idx) >= 0.1*size(fNames,1) % 10% of files
                    keepGenes = [keepGenes; U(ii)];
                end
            end
        end
        
        % convert back to table datatype and save to XLSX
        keepGenes = table(keepGenes);
        writetable(keepGenes,[srcDir 'finalGeneClusters_' num2str(eqClass) ...
            '.xlsx'],'WriteVariableNames',false,'WriteRowNames',false);
        
        eqClass = eqClass + 1; % increment equivalence class 1
        genes = []; % reinitialize to 0 so new vertcat can start
    end
end
fclose(fid);
