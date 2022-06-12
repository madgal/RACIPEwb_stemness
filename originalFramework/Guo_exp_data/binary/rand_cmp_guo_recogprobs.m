% calculate th p value
% the percent of predictions matched with experiments.

function rlt = rand_cmp_guo_recogprobs(num)

T = load('mouse_type.mat');
P = load('pred_type.mat');  % mean method
probs = [0.252484866601898,0.201503908527016,0.160189821388536,0.0542934010910993,0.0465959195874748,0.0458485912861520,0.0411030565727524,0.0359838577086914,0.0281542769598685,0.0222330169643524,0.0198041999850534,0.0189074060234661,0.0158807264031089,0.0136387414991406,0.00822061131455048];

% Rand
rlt  = zeros(1, num);
rlt0 = 0;
for j = 1:15
    if cmp_state(T.type, P.type(j,:))
        rlt0 = rlt0 + probs(j);
    end
end

for i = 1:num

    type = generate_type(15);
    
    for j = 1:15
        if cmp_state(T.type, type(j, :))
            rlt(i) = rlt(i) + probs(j);
%             display(j);
        end
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

function out = cmp_state(x, y)
out = 0;

for i = 1:size(x, 1)
    
    cnt = 1;
    for j = 1:size(x, 2)
        if (x(i, j) ~= 0 && y(j) ~= 0)
            if x(i, j) ~= y(j)
                cnt = 0;
            end
        end
    end
    
    if cnt == 1
        out = out + 1;
        break;
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