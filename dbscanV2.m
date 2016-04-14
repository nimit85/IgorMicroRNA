function [C, Core, Border, Outlier] = dbscanV2(D, radius, minpts)

% input data is in row format
Core = zeros(length(D),1); % vector stores core points
id = zeros(length(D),1); % vector stores cluster labels
Nbrs = cell(length(D),1); % cell array stores neighbors for each data point

for iter = 1 : length(D)
    % compute neighborhood of each x_i in D
    Nbrs{iter} = nbrFunc(iter, D, radius);
    if sum(Nbrs{iter}) >= minpts
        Core(iter) = 1;
    end
end

%%

k = 0; % cluster id

% run the DENSITYCONNECT algorithm to connect border points to the
% respective core points, and create final clusters
for iter = 1 : length(Core)
    if Core(iter)
        if ~id(iter)
            k = k + 1; % increase cluster label number
            id(iter) = k;
            id = densityConnect(iter, k, id, Nbrs, Core, D);
        end
    end
end

%%

% Collect all points within clusters
U = unique(id);
C = cell(numel(U),1);
for iter = 1 : numel(U)
    C{iter} = find(id == iter);
end

% Collect all outlier points (Outlier)
Outlier = find(~id);

% Collect all border points
setOnes = ones(size(D,1),1);
Border = (xor(setOnes, or(Core,~id)));

% plot the clusters and outliers
figure; hold on;
for iter = 1 : numel(C)
    thisClustPts = C{iter};
    if size(D,2) > 1 % multivariate data
        plot(D(thisClustPts,1),D(thisClustPts,2),'Color',rand(3,1)',...
            'Marker','o','LineStyle','none');
    else % univariate data
        plot(D(thisClustPts),'Color',rand(3,1)','Marker','o',...
            'LineStyle','none');
    end
end

if size(D,2) > 1 % multivariate data
    plot(D(Outlier,1),D(Outlier,2),'Color','k','Marker','*','LineStyle','none');
    % plot(D(Border,1),D(Border,2),'Color','m','Marker','+','LineStyle','none');
else % univariate data
    plot(D(Outlier),'Color','k','Marker','*','LineStyle','none');
    % plot(D(Border),'Color','m','Marker','+','LineStyle','none');
end
end

%%

% function calculates neighbors
function Nbrs = nbrFunc(thisPt, D, radius)
temp = D(thisPt,:);
if size(D,2) > 1
    tempVar = (repmat(temp,size(D,1),1) - D)';
    Nbrs = sqrt(sum(tempVar.^2))' < radius;
else
    Nbrs = abs(temp - D) < radius;
end
end
%%

% function combines core points and border points into same cluster, and
% identifies outliers
function id = densityConnect(thisIndex, k, id, Nbrs, Core, D)
thisPtNbrs = find(Nbrs{thisIndex} > 0);
for iter = 1 : length(thisPtNbrs)
    if ~id(thisPtNbrs(iter))
        id(thisPtNbrs(iter)) = k;
        if Core(thisPtNbrs(iter)) 
            id = densityConnect(thisPtNbrs(iter), k, id, Nbrs, Core, D);
        end
    end
end
end
