% calculate MSE fitting y = x
% the percent of predictions matched with experiments.
% without clustering before.
% shuffle the date to generate the random gene expression data
% want to know which predictions match which experiments.

function rlt = rand_cmp_guo_recogprobs_nocluster_ID

rng('shuffle');

mT = zeros(6, 9);
for i = 1:9
    T = importdata(['cluster_' num2str(i) '.dat']);
    mT(:, i) = mean(T)';
end

WD = importdata('../stem_solution_0_zscore.dat');

% wt 0.3532
rlt = zeros(size(WD.data, 1), 1);
for j = 1:size(WD.data, 1)
    WDv = WD.data(j, [2 3 5 6 7 9])';
    tmp = cmp_state(mT, WDv);
    if (tmp ~= 0)
        rlt(j) = tmp;
    end
end

end

function out = cmp_state(x, y)

% cut for each exp. gene state
% guassian
% cut = [0.4411, 0.4440, 0.1762, 0.1984, 0.2359, 0.8649, 0.2043, 0.2293, 0.2044]; %2.5%
% cut = [0.3995, 0.4009, 0.1617, 0.1805, 0.2151, 0.8317, 0.2017, 0.2105, 0.1828]; % 2.0%
% wt dist 1%
cut = [0.2479    0.2613    0.1523    0.1531    0.1841    0.5862    0.1559    0.1617    0.1600];

out = 0;

outx = [];
outy = [];

for i = 1:size(x, 2)
    
    tmp = sum(((x(:, i) - mean(x(:, i))) - (y - mean(y))).^2)/6;
    
    if tmp <= cut(i)
       outx = [outx i];
       outy = [outy abs(tmp - cut(i))/cut(i)];
    end
    
end

if ~isempty(outx)
    out = outx(outy == max(outy));
end

end