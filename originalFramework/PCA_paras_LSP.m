%% 9 genes
V = importdata('stem_solution_0_zscore.dat');

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data, 'Algorithm', 'svd');

% paras LSP
[N C]= hist3(SCORE(:, 1:2), [25 25]);
[x, y] = meshgrid(min(C{1}):0.1:max(C{1}), min(C{2}):0.1:max(C{2}));
Z = interp2(C{1}, C{2}, N'./sum(sum(N')), x, y);
surf(x, y, -log(Z));

%% 6 genes
V = importdata('stem_solution_0_zscore.dat');

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data(:, [2 3 5 6 7 9]), 'Algorithm', 'svd');

% paras LSP
[N C]= hist3(SCORE(:, 1:2), [25 25]);
[x, y] = meshgrid(min(C{1}):0.1:max(C{1}), min(C{2}):0.1:max(C{2}));
Z = interp2(C{1}, C{2}, N'./sum(sum(N')), x, y);
surf(x, y, -log(Z));

G = importdata('Guo_data_all_extracted_data.dat');
[COEFF2, SCORE2, LATENT2, TSQUARED2, EXPLAINED2] = pca(G.data(:, 2:end), 'Algorithm', 'svd');

[N C]= hist3(SCORE2(:, 1:2), [10 10]);
[x, y] = meshgrid(min(C{1}):0.1:max(C{1}), min(C{2}):0.1:max(C{2}));
Z = interp2(C{1}, C{2}, N'./sum(sum(N')), x, y);
figure;
surf(x, y, -log(Z));

figure;
plot(SCORE2(:, 1), SCORE2(:, 2), '.');
axis square

wtx = zeros(size(G.data, 1), 1);
wty = zeros(size(G.data, 1), 1);

for i = 1:6
    wtx = wtx + COEFF(i, 1).*G.data(:, i+1);
    wty = wty + COEFF(i, 2).*G.data(:, i+1);
end

[N C]= hist3([wtx, wty], [10 10]);
[x, y] = meshgrid(min(C{1}):0.1:max(C{1}), min(C{2}):0.1:max(C{2}));
Z = interp2(C{1}, C{2}, N'./sum(sum(N')), x, y);
figure;
surf(x, y, -log(Z));
