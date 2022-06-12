% calculate th p value for recognized number of experimental gene states by
% % calculate MSE fitting y = x

function rlt = rand_cmp_guo_num(num)

rng('shuffle');

mT = zeros(6, 9);
for i = 1:9
    T = importdata(['cluster_' num2str(i) '.dat']);
    mT(:, i) = mean(T)';
end

W = importdata('../Wt_average.dat');
WD = importdata('../stem_solution_0_zscore.dat');
W2 = W.data(:, [2 3 5 6 7 9])';

% cut for each exp. gene state
% guassian
% cut = [0.4411, 0.4440, 0.1762, 0.1984, 0.2359, 0.8649, 0.2043, 0.2293, 0.2044]; %2.5%
% cut = [0.3995, 0.4009, 0.1617, 0.1805, 0.2151, 0.8317, 0.2017, 0.2105, 0.1828]; % 2.0%
% wt dist
cut = [0.2479    0.2613    0.1523    0.1531    0.1841    0.5862    0.1559    0.1617    0.1600];

% Rand
rlt  = zeros(1, num);
rlt0 = 0;
for j = 1:9
    rlt0 = rlt0 + cmp_state(W2, mT(:, j), cut(j));
end

pos  = [2 3 5 6 7 9];
for i = 1:num

    % type = randn([6, 15]);
    
    type = zeros(6, 15);
    
    for j = 1:6
        type(j, :) = WD.data(randi(size(WD.data, 1), [15, 1]), pos(j))';
    end
    
%     clc;
    
    for j = 1:9
        rlt(i) = rlt(i) + cmp_state(type, mT(:, j), cut(j));
    end

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

function out = cmp_state(x, y, cut)

out = 0;

for i = 1:size(x, 2)
    
    tmp = sum(((x(:, i) - mean(x(:, i))) - (y - mean(y))).^2)/6;
    
    if tmp <= cut
        out = 1;
%         display(tmp);
%         display(x(:, i));
%         display(y);
    end
    
end



end

function type = generate_type(num2)

ngene = 6;

type = zeros(num2, ngene);

cnt = 0;

while (cnt < num2)
    if cnt == 0
        type(1,:) = randn([1 6]);
        cnt = cnt + 1;
    else
        tmp = randn([1 6]);
        cnt2 = 0;
        for i = 1:cnt
           tmpcorr = corr(type(i, :)', tmp');
           if tmpcorr < 0.72
               cnt2 = 1;
           end
        end
        
        if cnt2 == 0
            cnt = cnt + 1;
            type(cnt, :) = tmp;
        end
    end

end

end