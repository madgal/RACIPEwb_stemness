function hist_data_guo

V0 = importdata('Guo_data.dat');
figure;

% gene_state_hist_cmp_2.fig
for num = 1:9
    V  = importdata(['cluster_' num2str(num) '.dat']);

    cnt = 50;
    
    x = linspace(-3, 3, cnt + 1);

    for i = 1:6
        probs0 = zeros(1, cnt);
        probs  = zeros(1, cnt);

        for j = 1:cnt

            probs0(j) = length(find(V0(:,i) >= x(j) & V0(:,i) < x(j+1)))/size(V0, 1);
            probs (j) = length(find(V (:,i) >= x(j) & V (:,i) < x(j+1)))/size(V0, 1);

        end

        subplot(12, 6, (num - 1)*6 + i);
        hold on;
        bar(x(1:cnt) + diff(x)/2, probs0./max(probs0), 'FaceColor', [0.8960 0.8960 0.8960], 'EdgeColor', [0.8960 0.8960 0.8960]);
        bar(x(1:cnt) + diff(x)/2, probs ./max(abs(probs)) , 'FaceColor', [1      0      0     ], 'EdgeColor', [1      0      0     ]);
        ylim([-1 1]);
        hold off;

    end
end


% gene_state_hist_cmp.fig
% for num = 1:12
%     V  = importdata(['Wt_cluster_' num2str(num) '.dat']);
% 
%     cnt = 50;
%     
%     x = linspace(-3, 3, cnt + 1);
% 
%     for i = 1:9
%         probs0 = zeros(1, cnt);
%         probs  = zeros(1, cnt);
% 
%         for j = 1:cnt
% 
%             probs0(j) = length(find(V0(:,i) >= x(j) & V0(:,i) < x(j+1)))/size(V0, 1);
%             probs (j) = -1.*length(find(V (:,i) >= x(j) & V (:,i) < x(j+1)))/size(V, 1);
% 
%         end
% 
%         subplot(12, 9, (num - 1)*9 + i);
%         hold on;
%         bar(x(1:cnt) + diff(x)/2, probs0./max(probs0), 'FaceColor', [0.8960 0.8960 0.8960], 'EdgeColor', [0.8960 0.8960 0.8960]);
%         bar(x(1:cnt) + diff(x)/2, probs ./max(abs(probs)) , 'FaceColor', [1      0      0     ], 'EdgeColor', [1      0      0     ]);
%         ylim([-1 1]);
%         hold off;
% 
%     end
% end