V = importdata('stem_solution_0_zscore.dat');

% up layer
fid1 = fopen('up_data.dat', 'w');
fprintf(fid1, '%s\t%s\t%s\t%s\t%s\n', 'Name', 'Cdx2', 'Oct4', 'Oct4-Sox2', 'Sox2');
for i = 1:length(V.data)
  fprintf(fid1, '%s\t%f\t%f\t%f\t%f\n', V.textdata{i+1, 1}, V.data(i, [2, 7, 8, 9]));
end
fclose(fid1);

% bottom layer
fid2 = fopen('bottom_data.dat', 'w');
fprintf(fid2, '%s\t%s\t%s\t%s\t%s\t%s\n', 'Name', 'Gcnf', 'Gata6', 'Pbx1', 'Klf4', 'Nanog');
for i = 1:length(V.data)
  fprintf(fid2, '%s\t%f\t%f\t%f\t%f\t%f\n', V.textdata{i+1, 1}, V.data(i, [1, 3, 4, 5, 6]));
end
fclose(fid2);