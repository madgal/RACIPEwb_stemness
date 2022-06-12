V = importdata('stem_solution_0_zscore.dat');
% V = importdata('stem_solution_0_raw.dat');
ID = importdata('stem_solution_0_result_wt_avg.dat');

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data, 'Algorithm', 'svd');

figure;
plot(cumsum(EXPLAINED), 'o--');
xlabel('Principle components');
ylabel('Variance explained');


% PC1 and PC2
figure('Color', [1 1 1]);
% subplot(2, 3, 1);
hist3(SCORE(:, 1:2), [100 100]);
set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
axis square;

%{
subplot(2, 3, 2);
cmap = jet(15);
hold on;
for j = 1:15
    plot(SCORE(ID.data == j, 1), SCORE(ID.data == j, 2), '.', 'MarkerSize', 2, 'Color', cmap(j, :));  
end

for j = 1:15
    text(mean(SCORE(ID.data == j, 1)), mean(SCORE(ID.data == j, 2)), num2str(j), 'FontSize', 15);
end

colormap(cmap);
hold off;
axis square;

subplot(2, 3, 3);
hold on;
% plot(COEFF(:, 1), COEFF(:, 2), 'ko', 'MarkerSize', 8);
text(COEFF(:, 1), COEFF(:, 2),  cellstr(['Gcnf     '; 'Cdx2     '; 'Gata6    '; 'Pbx1     '; 'Klf4     '; 'Nanog    '; 'Oct4     '; 'Oct4-Sox2'; 'Sox2     ']));
plot([0 0], [-0.4 0.8], 'r--', 'LineWidth', 1);
plot([-0.4 0.4], [0 0], 'r--', 'LineWidth', 1);
xlim([-0.4, 0.4]);
ylim([-0.4, 0.8]);
hold off;
axis square;

% PC3 and PC4
subplot(2, 3, 4);
hist3(SCORE(:, 3:4), [100 100]);
set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
axis square;

subplot(2, 3, 5);
cmap = jet(15);

hold on;
for j = 1:15
    plot(SCORE(ID.data == j, 3), SCORE(ID.data == j, 4), '.', 'MarkerSize', 2, 'Color', cmap(j, :));
end

for j = 1:15
    text(mean(SCORE(ID.data == j, 3)), mean(SCORE(ID.data == j, 4)), num2str(j), 'FontSize', 15);
end

colormap(cmap);
hold off;
axis square;

subplot(2, 3, 6);
hold on;
% plot(COEFF(:, 3), COEFF(:, 4), 'ko', 'MarkerSize', 8);
text(COEFF(:, 3), COEFF(:, 4), cellstr(['Gcnf     '; 'Cdx2     '; 'Gata6    '; 'Pbx1     '; 'Klf4     '; 'Nanog    '; 'Oct4     '; 'Oct4-Sox2'; 'Sox2     ']));
plot([0 0], [-0.4 0.8], 'r--', 'LineWidth', 1);
plot([-0.4 0.4], [0 0], 'r--', 'LineWidth', 1);
xlim([-0.4, 0.4]);
ylim([-0.4, 0.8]);
hold off;
axis square;
%}