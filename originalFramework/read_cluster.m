% use ./colorByThreshold.pl -t 0.8 -c ID.color stem_solution_0_zscore.gtr
% then use this to see the probability distributution. 
% remove first line of input file.

% thrd = 0.8, 53 clusters

function cnt = read_cluster(num)

fid = fopen('stem_solution_0_zscore_gtr.dat','r');

frewind(fid);

data = textscan(fid, '%s\t%s\t%s\t%f\t%d\n');

cnt = zeros(1,num);

for i = 1:num
    idx = find(data{5} == i);
    tmpcnt = 0;
    for j = 1:length(idx)
        for h = 2:3
            if strncmpi(data{h}(idx(j)), 'GENE', 4)
                tmpcnt = tmpcnt + 1;
            end
        end
    end
    cnt(i) = tmpcnt;
end

fclose(fid);

figure; 
bar(sort(cnt./sum(cnt), 2, 'descend'));
hold
plot([0 num+1], [0.005 0.005], 'r');

end

