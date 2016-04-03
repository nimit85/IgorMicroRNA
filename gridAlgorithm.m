load('gene2.mat', 'genes');

X = zscore(genes);
Q=X;

% grid method
l = -2;
r = 18;
b = -5;
t = 25;
g = 0.6;

figure
plot(Q(:,1), Q(:, 2), 'b.');
hold all
numRows = floor((r-l)/g);
numCols = floor((t-b)/g);
% fill in the grid
Grid = zeros(numRows, numCols);
for i = 1:numRows
    for j = 1:numCols
        gl = (i-1)*g+l;
        gb = (j-1)*g+b;
        gr = (i)*g+l;
        gt = (j)*g+b;
        
        glr = Q(:,1) > gl & Q(:,1) <= gr;
        gbt = Q(:,2) > gb & Q(:,2) <= gt;
        
        ing = glr & gbt;
        
        plot([gl gr gl gr], [gb gt gt gb], 'k+');
        plot(Q(ing,1), Q(ing, 2), 'color', rand(1,3), ...
            'marker', '.', 'linestyle', 'none');
        if sum(glr & gbt) > 1
            Grid(i,j) = sum(glr & gbt)/g^2;
        end
    end
end
figure
image(imrotate(Grid, 90));

% function to measure distance between grid cells
dsim = @(x, y) abs(x-y)/(x+y);
thresh = 0.95;

% start with every cell in a different cluster
Class = reshape(randperm(length(Grid(:))), size(Grid));
Class(Grid==0) = 0;

% color code the clusters
Grid2 = rand(numRows, numCols, 3);
for i = 1:numRows
    for j = 1:numCols
        if Grid(i,j) == 0
            Grid2(i,j,:) = [1,1,1]';
        end
    end
end

figure
for t = 1:50
% merge similar cells (take from low,left and pass upper right
for i = 2:1:numRows
    for j = 2:1:numCols
%         lower left
        if dsim(Grid(i,j),Grid(i-1, j-1)) <= thresh
            Class(i,j) = Class(i-1,j-1);
            Grid2(i,j,:) = Grid2(i-1,j-1,:);
        end
%         lower
        if dsim(Grid(i,j),Grid(i, j-1)) <= thresh
            Class(i,j-1) = Class(i,j);
            Grid2(i,j-1,:) = Grid2(i,j,:);
        end
%         left
        if dsim(Grid(i,j),Grid(i-1, j)) <= thresh
            Class(i-1,j) = Class(i,j);
            Grid2(i-1,j,:) = Grid2(i,j,:);
        end        
    end
end
image(imrotate(Grid2, 90));
title(['Iteration: ' num2str(t)])
pause(0.1)
end





