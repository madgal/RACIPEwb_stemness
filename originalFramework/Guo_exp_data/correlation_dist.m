% 0.72 is a good cutoff for pearson correlation (5%)
% MSE for fitting y = x, where the data should be centered before fitting.
% cut off is 2.5%
function [cut, coeffw] = correlation_dist(num, nrun)

rng('shuffle');

T = importdata(['cluster_' num2str(num) '.dat']);

mT = mean(T)';

coeff = zeros(1, nrun);

WD = importdata('../stem_solution_0_zscore.dat');

for i = 1:nrun
    
%     R = randn([6, 1]);
    R = zeros(6, 1);
    pos  = [2 3 5 6 7 9];
    for j = 1:6
        R(j) = WD.data(randi(size(WD.data, 1), 1), pos(j));
    end
    
    X = (mT - mean(mT));
    Y = (R - mean(R));
    
%     coeff(i) = corr(mT, R); % pearson correlation
%     coeff(i) = sum((X - Y).^2)/6;  % MSE
    
    coeff(i) = 1 - sum((Y - X).^2)/sum((Y - mean(Y)).^2); % RSS
    
end

W = importdata('../Wt_average.dat');
W2 = W.data(:, [2 3 5 6 7 9])';

coeffw = zeros(1, size(W2, 2));

for i = 1:size(W2, 2)
   
   X = (mT - mean(mT));
   Y = (W2(:, i) - mean(W2(:, i)));
   %coeffw(i) = corr(mT, W2(:, i)); % pearson correlation
%    coeffw(i) = sum((X - Y).^2)/6;  % MSE
   coeffw(i) = 1 - sum((Y - X).^2)/sum((Y - mean(Y)).^2); % RSS
end
    
scoeff = sort(coeff, 'ascend');
p = 0.01;
cut = scoeff(nrun*p);

% plot
%%{
% figure('Color', [1 1 1]);
hold on;
histogram(coeff, 'Normalization', 'cdf', 'BinWidth', 0.05, 'FaceColor', [0.6000 0.8 1], 'EdgeColor', [0.6000 0.8 1]);
for i = 1:size(W2, 2)
    plot([coeffw(i) coeffw(i)], [0 1], 'k-', 'LineWidth', 3);
    if (coeffw(i) <= 2 && coeffw(i) >= 0) % need to change boundary
        text(coeffw(i), 1.1, num2str(i));
    end
end
ylim([0 1.2]);

%MSE
% plot([0 2], [p p], 'r--'); %MSE
% plot([scoeff(nrun*p) scoeff(nrun*p)], [0 1], 'r-');
% xlim([0 2]);

% RSS
plot([0 1], [1-p 1-p], 'r--');
% plot([scoeff(nrun*(1-p)) scoeff(nrun*(1-p))], [0 1], 'r-');
xlim([0 1]);
%}
