function separate_clusters

% get numclt by cnt = read_cluster(38); [Y,I] = sort(cnt,2,'descend')
numclt = [2,53,51,4,39,5,31,36,17,21,19,6,37,40,28];

gtr = importdata('stem_solution_0_zscore_gtr.dat');
cdt = importdata('stem_solution_0_zscore_cdt.dat');

for i = 1:length(numclt)
    
   idx = find(gtr.data(:,2) == numclt(i));
   fid = fopen(['Wt_cluster_' num2str(i) '.dat'], 'w');
   
   for j = 1:length(idx)
      
       for h = 1:2
           if strncmp(gtr.textdata{idx(j), h+1}, 'GENE', 4)
               
               for m = 1:length(cdt.data)
                   if strcmp(cdt.textdata{m,1}, gtr.textdata{idx(j), h+1})
                       break;
                   end
               end
               
               fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', cdt.data(m, 2:end));
           end
       end
       
   end
   
   fclose(fid);
end

end