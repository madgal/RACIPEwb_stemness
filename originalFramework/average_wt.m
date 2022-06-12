function average_wt(num)

%% wt data
cmpdata = cell(1, num);
wt_data = zeros(num,9);

for i = 1:num
    cmpdata{i}  = importdata(['Wt_cluster_' num2str(i) '.dat']);
    wt_data(i, :) = mean(cmpdata{i});
end

fid = fopen('Wt_average.dat', 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'name', 'Gcnf', 'Cdx2', 'Gata6', 'Pbx1', 'Klf4', 'Nanog', 'Oct4', 'Oct4-Sox2', 'Sox2');

for i = 1:num
    fprintf(fid, '%s\t', ['T_' num2str(i)]);
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', wt_data(i, :));
end

fclose(fid);
