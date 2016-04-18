function [C, Core, Border, Noise] = dbscanV2(D, radius, minpts)

D = zscore(D);
% input data is in row format
Core = zeros(length(D),1);
id = zeros(length(D),1);
Nbrs = cell(length(D),1);
for iter = 1 : length(D)
    % compute neighborhood of each x_i in D
    Nbrs{iter} = nbrFunc(iter, D, radius);
    if sum(Nbrs{iter}) >= minpts
        Core(iter) = 1;
    end
end

disp('Core Points calculated');

%%

k = 0; % cluster id

for iter = 1 : length(Core)
    if Core(iter)
        if ~id(iter)
            k = k + 1;
            id(iter) = k;
            id = densityConnect(iter, k, id, Nbrs);
        end
    end
end

disp('Clusters calculated');
%%

% Collect all points within clusters
U = unique(id);
C = cell(numel(U),1);

minVal = Inf;
for iter = 1 : numel(U)
    C{iter} = find(id == iter);
    if length(C{iter}) < minVal && ~isempty(C{iter})
        minVal = length(C{iter});
    end
end

% Collect all outlier points (noise)
Noise = find(~id);

% Collect all border points
setOnes = ones(size(D,1),1);
Border = (xor(setOnes, or(Core,~id)));

% random sampling cntrOnS points from each cluster - samplesClust contains
% the randomly sampled points
samplesClust = zeros(length(U),minVal);
for j = 1 : size(C,1)
    if ~isempty(C{j})
        F = C{j};
        onlyFew = F(randperm(size(F,1)));
        samplesClust(j,:) = onlyFew(1:minVal);
    end
end
samplesClust(end,:) = [];

load('geneNames.mat');
C = B1H_0I_M0_R1_Jcel(samplesClust);
T = cell2table(C);
writetable(T,'DBScanOutliers.txt','WriteRowNames',false,'WriteVariableNames',false);

%%
% plot the clusters
figure; hold on;
for iter = 1 : numel(C)
    thisClustPts = C{iter};
    if size(D,2) > 1
        plot(D(thisClustPts,1),D(thisClustPts,2),'Color',rand(3,1)',...
            'Marker','o','LineStyle','none');
    else
        plot(D(thisClustPts),'Color',rand(3,1)','Marker','o',...
            'LineStyle','none');
    end
end
if size(D,2) > 1
    plot(D(Noise,1),D(Noise,2),'Color','k','Marker','+','LineStyle','none');
    plot(D(Border,1),D(Border,2),'Color','m','Marker','d','LineStyle','none');
else
    plot(D(Noise),'Color','k','Marker','+','LineStyle','none');
    plot(D(Border),'Color','m','Marker','d','LineStyle','none');
end
legend('show');

end

function Nbrs = nbrFunc(thisPt, D, radius)
temp = D(thisPt,:);
if size(D,2) > 1
    Nbrs = (pdist2(temp,D) < radius);
else
    Nbrs = abs(temp - D) < radius;
end
end
%%

function id = densityConnect(thisIndex, k, id, Nbrs)
thisPtNbrs = find(Nbrs{thisIndex} > 0);
for iter = 1 : length(thisPtNbrs)
    if ~id(thisPtNbrs(iter))
        id(thisPtNbrs(iter)) = k;
%         if Core(thisPtNbrs(iter)) 
%             id = densityConnect(thisPtNbrs(iter), k, id, Nbrs, Core, D);
%         end
    end
end
end