% calculate MSE fitting y = x
% the percent of predictions matched with experiments.
% without clustering before.
% shuffle the date to generate the random gene expression data

function rlt = rand_cmp_guo_recogprobs_nocluster(num)

rng('shuffle');

mT = zeros(6, 9);
for i = 1:9
    T = importdata(['cluster_' num2str(i) '.dat']);
    mT(:, i) = mean(T)';
end

% W = importdata('../Wt_average.dat');
WD = importdata('../stem_solution_0_zscore.dat');
% W2 = W.data(:, [2 3 5 6 7 9])';

% probs = [0.252484866601898,0.201503908527016,0.160189821388536,0.0542934010910993,0.0465959195874748,0.0458485912861520,0.0411030565727524,0.0359838577086914,0.0281542769598685,0.0222330169643524,0.0198041999850534,0.0189074060234661,0.0158807264031089,0.0136387414991406,0.00822061131455048];

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
%%{
pos  = [2 3 5 6 7 9];
for i = 1:num

%     type = randn([6, 15]);
    type = zeros(6, size(WD.data, 1));
    
    for j = 1:6
        type(j, :) = WD.data(randi(size(WD.data, 1), [size(WD.data, 1), 1]), pos(j))';
    end

    for j = 1:size(WD.data, 1)
        if cmp_state(mT, type(:, j))
            rlt(i) = rlt(i) + 1;
%             display(j);
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
%}

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

function type = generate_type(num2)

ngene = 9;

type = zeros(num2, ngene);

cnt = 0;

while (cnt < num2)
    if cnt == 0
        type(1,:) = randi([0 1], [1, ngene]);
        cnt = cnt + 1;
    else
        tmp = randi([0 1], [1, ngene]);
        cnt2 = 0;
        for i = 1:cnt
           if isequal(type(i, :), tmp)
               cnt2 = 1;
           end
        end
        
        if cnt2 == 0
            cnt = cnt + 1;
            type(cnt, :) = tmp;
        end
    end

end

for i = 1:ngene
    for j = 1:num2
        if type(j, i) == 0
            type(j, i) = -1;
        end
    end
end

% ramdom assgn 0
% rw = randi([1 num2], [1, 2]);
% cn = randi([1 9], [1, 2]); 
% for i = 1:2
%     type(rw(i), cn(i)) = 0;
% end

type = type(:,1:6);

end