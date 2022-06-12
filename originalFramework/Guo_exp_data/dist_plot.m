function dist_plot

figure;

for i = 1:9
    
    V = importdata(['cluster_' num2str(i) '.dat']);
    
    for j = 1:6
        subplot(9, 6, (i-1)*6 + j);
        histogram(V(:,j), [-4 0 4], 'Normalization', 'Probability');
        xlim([-3 3]);
    end
end

