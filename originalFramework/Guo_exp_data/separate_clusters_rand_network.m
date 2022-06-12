% for .cdt, need to remove the first 2 lines.
% for .gtr, need to remove the first line.

function separate_clusters_rand_network

numclt = 1:1:12;

gtr = importdata('Guo_data_all_extracted_cnt.dat');
cdt = importdata('Guo_data_all_extracted_data.dat');
cnt = 0;

for i = 1:length(numclt)
    
   idx = find(gtr.data(:,2) == numclt(i));
   if length(idx)/length(cdt.data) >= 0.01 
       cnt = cnt + 1;
       fid = fopen(['cluster_' num2str(cnt) '.dat'], 'w');

       for j = 1:length(idx)
           for h = 1:2
               if strncmp(gtr.textdata{idx(j), h+1}, 'GENE', 4)

                   for m = 1:length(cdt.data)
                       if strcmp(cdt.textdata{m,1}, gtr.textdata{idx(j), h+1})
                           break;
                       end
                   end

                   fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\n', cdt.data(m, 2:end));
               end
           end

       end
       fclose(fid);

   end
end

end