V = importdata('stem_solution_0_result_wt_avg.dat');

stb = zeros(3, 6);

figure('Color', [1 1 1]);

for i = 1:3
    idx = find(V.data == i);
    
    for j = 1:length(idx)
        tmp = strsplit(V.textdata{idx(j)}, '_');
        num = str2num(tmp{1});
        stb(i, num) = stb(i, num) + 1;
    end
   
    subplot(1, 3, i);
    
    bar(stb(i, :)./sum(stb(i, :)));
end