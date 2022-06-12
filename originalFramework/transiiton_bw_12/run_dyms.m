clc;clear;

Edge{1} = linspace(-5, 5, 100);
Edge{2} = linspace(-3, 3, 100);

V = importdata('../stem_solution_0_zscore.dat');

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data, 'Algorithm', 'svd');

friv   = importdata('RIV_12.dat');
idx = find(friv(:, 3) == 1);

y = cell(1, 2);

for i = 1:2
    y{i} = zeros(1501, length(idx));
end

for i = 1:length(idx)
    ytmp = dyms(1500, idx(i));
    
    G = zeros(size(ytmp, 1), 2);
    pos = [2 3 1 6 4 5 7 8 9];

    for j = 1:9
        G(:, 1) = G(:, 1) + COEFF(j, 1).*ytmp(:, pos(j));
        G(:, 2) = G(:, 2) + COEFF(j, 2).*ytmp(:, pos(j));
    end
    
    for j = 1:2
        y{j}(:, i) = G(1:10:end, j);
    end
end

figure;
hold on;
hist3(SCORE(:, 1:2), 'Edges',Edge);
set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
axis square;
xlim([-5 5]);
ylim([-3 3]);
view(2);

for i = 1:length(idx)
    plot3(y{1}(:, i), y{2}(:, i), ones(1501, 1).*1000, 'w.-');
end

hold off;

% plot
x1 = []; 
for i = 1:length(idx); 
    x1 = [x1; y{1}(:, i)]; 
end;
x2 = []; 
for i = 1:length(idx); 
    x2 = [x2; y{2}(:, i)]; 
end;
x = [x1, x2];
Edge{1} = linspace(-5, 5, 100);
Edge{2} = linspace(-3, 3, 100);
figure;
hist3(x, 'Edges',Edge);
set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
axis square;


Edge{1} = linspace(-5, 5, 100);
Edge{2} = linspace(-3, 3, 100);
figure;
hold on;
hist3(SCORE(:, 1:2), 'Edges',Edge);
set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
axis square;
xlim([-5 5]);
ylim([-3 3]);
view(2);

num = 100;
Edge{1} = linspace(-5, 5, num);
Edge{2} = linspace(-3, 3, 100);
ym = zeros(num-1, 1);
for i = 1:num-1
   tmp = 0;
   cnt = 0;
   for j = 1:size(y{1}, 2) 
       idx = find(y{1}(:, j) <= Edge{1}(i+1) & y{1}(:, j) > Edge{1}(i));
       cnt = cnt + length(idx);
       tmp = tmp + sum(y{2}(idx, j));
   end
   ym(i) = tmp/cnt;
end

plot3(Edge{1}(1:num-1) + diff(Edge{1})/2, ym, ones(num-1,1).*100, 'wo--');


