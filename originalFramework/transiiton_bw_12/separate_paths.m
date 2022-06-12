clc;clear;

load('transition_12_a_nodecrease.mat');

IDS = importdata('../stem_solution_0_result_wt_min.dat');

ID = size(1, length(idx));
cnt = zeros(1, 3);
for i = 1:length(idx)
    if length(find(y{1}(:, i) <= -1.8)) >= 10
        if ~isempty(find(y{2}(:, i) >= 1.5))
            ID(i) = 1;
            cnt(1) = cnt(1) + 1;
        elseif ~isempty(find(y{2}(:, i) >= 0.5))
            ID(i) = 2;
            cnt(2) = cnt(2) + 1;
        else
            ID(i) = 3;
            cnt(3) = cnt(3) + 1;
        end
    end
end

% separate plot of trajs
figure('Color', [1 1 1]);
for i = 1:3
    subplot(1, 3, i);
%     hist3(SCORE(:, 1:2), 'Edges',Edge);
%     set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
%     axis square;
%     xlim([-5 5]);
%     ylim([-3 3]);
%     view(2);

    cmap = jet(15);
    hold on;
    for j = 1:15
        plot(SCORE(IDS.data == j, 1), SCORE(IDS.data == j, 2), '.', 'MarkerSize', 2, 'Color', cmap(j, :));  
    end

    for j = 1:15
        text(mean(SCORE(IDS.data == j, 1)), mean(SCORE(IDS.data == j, 2)), num2str(j), 'FontSize', 15);
    end

    colormap(cmap);
    xlim([-5 5]);
    ylim([-3 3]);
    view(2);
    axis square;
    
    ID2 = find(ID == i);
%     for j = 1:length(ID2)
%         plot3(y{1}(:, ID2(j)), y{2}(:, ID2(j)), ones(size(y{2}(:, ID2(j)))).*1000, 'w.-');
%     end
%     
    text(-4, -2, num2str(cnt(i)./sum(cnt)), 'Color', [1 1 1]);
    
    display(num2str(cnt(i)./sum(cnt)));
    
    % average path
    num = 100;
    Edge{1} = linspace(-5, 5, num);
    Edge{2} = linspace(-3, 3, 100);
    ym = zeros(num-1, 1);
    for h = 1:num-1
       tmp = 0;
       cnt2 = 0;
       for j = 1:length(ID2)
           idx = find(y{1}(:, ID2(j)) <= Edge{1}(h+1) & y{1}(:, ID2(j)) > Edge{1}(h));
           cnt2 = cnt2 + length(idx);
           tmp = tmp + sum(y{2}(idx, ID2(j)));
       end
       ym(h) = tmp/cnt2;
    end

    plot3(Edge{1}(1:num-1) + diff(Edge{1})/2, ym, ones(num-1,1).*1300, 'Ko--');
    
    hold off;
end



% % plot all individual trajs.
% figure;
% hold on;
% hist3(SCORE(:, 1:2), 'Edges',Edge);
% set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
% axis square;
% xlim([-5 5]);
% ylim([-3 3]);
% view(2);
% 
% for i = 1:length(idx)
%     plot3(y{1}(:, i), y{2}(:, i), ones(1501, 1).*1000, 'w.-');
% end
% 
% hold off;

% % heatmap plot
% x1 = []; 
% for i = 1:length(idx); 
%     x1 = [x1; y{1}(:, i)]; 
% end;
% x2 = []; 
% for i = 1:length(idx); 
%     x2 = [x2; y{2}(:, i)]; 
% end;
% x = [x1, x2];
% Edge{1} = linspace(-5, 5, 100);
% Edge{2} = linspace(-3, 3, 100);
% figure;
% hist3(x, 'Edges',Edge);
% set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
% axis square;
% 
% 
% Edge{1} = linspace(-5, 5, 100);
% Edge{2} = linspace(-3, 3, 100);
% figure;
% hold on;
% hist3(SCORE(:, 1:2), 'Edges',Edge);
% set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
% axis square;
% xlim([-5 5]);
% ylim([-3 3]);
% view(2);
% 
% % average path
% num = 100;
% Edge{1} = linspace(-5, 5, num);
% Edge{2} = linspace(-3, 3, 100);
% ym = zeros(num-1, 1);
% for i = 1:num-1
%    tmp = 0;
%    cnt = 0;
%    for j = 1:size(y{1}, 2) 
%        idx = find(y{1}(:, j) <= Edge{1}(i+1) & y{1}(:, j) > Edge{1}(i));
%        cnt = cnt + length(idx);
%        tmp = tmp + sum(y{2}(idx, j));
%    end
%    ym(i) = tmp/cnt;
% end
% 
% plot3(Edge{1}(1:num-1) + diff(Edge{1})/2, ym, ones(num-1,1).*100, 'wo--');