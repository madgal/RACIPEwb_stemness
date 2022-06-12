% calculate MSE fitting y = x
% the percent of predictions matched with experiments.
% without clustering before.
% shuffle the date to generate the random gene expression data
% but assume that there are some existing clusters internally.

function rlt = rand_cmp_guo_recogprobs_nocluster_2(num)

rng('shuffle');

mT = zeros(6, 9);
for i = 1:9
    T = importdata(['cluster_' num2str(i) '.dat']);
    mT(:, i) = mean(T)';
end

WD = importdata('../stem_solution_0_zscore.dat');

% wt 0.3532
rlt  = zeros(1, num);
rlt0 = 0;
for j = 1:size(WD.data, 1)
    WDv = WD.data(j, [2 3 5 6 7 9])';
    if cmp_state(mT, WDv)
        rlt0 = rlt0 + 1;
    end
end
rlt0 = rlt0/size(WD.data, 1);

%rand

for i = 1:num

    type = generate_type(15, WD.data);

    for j = 1:size(WD.data, 1)
        if cmp_state(mT, type(:, j))
            rlt(i) = rlt(i) + 1;
        end
    end
    
    rlt(i) = rlt(i)/size(WD.data, 1);

end

figure;
subplot(1,2,1);
hold on;
histogram(rlt, 'Normalization', 'PDF');
plot([rlt0, rlt0], [0 1], 'r-', 'LineWidth', 3);
hold off;
subplot(1,2,2);
hold on;
histogram(rlt, 'Normalization', 'CDF');
plot([rlt0, rlt0], [0 1], 'r-', 'LineWidth', 3);
hold off;

end

function out = cmp_state(x, y)

% cut for each exp. gene state
% guassian
% cut = [0.4411, 0.4440, 0.1762, 0.1984, 0.2359, 0.8649, 0.2043, 0.2293, 0.2044]; %2.5%
% cut = [0.3995, 0.4009, 0.1617, 0.1805, 0.2151, 0.8317, 0.2017, 0.2105, 0.1828]; % 2.0%
% wt dist 1%
cut = [0.2479    0.2613    0.1523    0.1531    0.1841    0.5862    0.1559    0.1617    0.1600];

out = 0;

for i = 1:size(x, 2)
    
    tmp = sum(((x(:, i) - mean(x(:, i))) - (y - mean(y))).^2)/6;
    
    if tmp <= cut(i)
        out = 1;
    end
    
end

end

function type = generate_type(num, WD)

    type = [];
    V = importdata('../stem_solution_0_result_wt_min.dat');
    
    for i = 1:num
       
       data = WD(V.data == i, [2 3 5 6 7 9]);
       
       ci = randperm(6, 6);
       
       % only shuffle columns
       typetmp = data(:, ci);
       
       % shuffle rows
%        typetmp = zeros(size(data));
%        for j = 1:6
%           ri = randi(size(data, 1), [size(data, 1) 1]);
%           typetmp(:, j) =  data(ri, ci(j));
%        end
        
       type = [type; typetmp];
       
    end
    
    type = type';

end