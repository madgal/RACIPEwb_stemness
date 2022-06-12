function classify_states_wt_min(casenum)

name = 'stem_solution';

source = importdata([name '_' num2str(casenum) '_zscore.dat']);

cmpdata{1}  = importdata('Wt_cluster_1.dat');
cmpdata{2}  = importdata('Wt_cluster_2.dat');
cmpdata{3}  = importdata('Wt_cluster_3.dat');
cmpdata{4}  = importdata('Wt_cluster_4.dat');
cmpdata{5}  = importdata('Wt_cluster_5.dat');
cmpdata{6}  = importdata('Wt_cluster_6.dat');
cmpdata{7}  = importdata('Wt_cluster_7.dat');
cmpdata{8}  = importdata('Wt_cluster_8.dat');
cmpdata{9}  = importdata('Wt_cluster_9.dat');
cmpdata{10} = importdata('Wt_cluster_10.dat');
cmpdata{11} = importdata('Wt_cluster_11.dat');
cmpdata{12} = importdata('Wt_cluster_12.dat');
cmpdata{13} = importdata('Wt_cluster_13.dat');
cmpdata{14} = importdata('Wt_cluster_14.dat');
cmpdata{15} = importdata('Wt_cluster_15.dat');

fid = fopen([name '_' num2str(casenum) '_result_wt_min.dat'], 'w');

for i = 1:15
    for j = 1:9
        for h = 1:length(cmpdata{i}(:,j))
            if cmpdata{i}(h,j) > 1
                cmpdata{i}(h,j) = 1;
            elseif cmpdata{i}(h,j) < -1
                cmpdata{i}(h,j) = -1;
            end
        end
    end
end

for i = 1:length(source.data)
   tmpM   = zeros(1,15);
   
   for j = 1:15
       [tmpM(j)] = dsnt(source.data(i,:), cmpdata{j});  
   end
   
   idx = find(tmpM == min(tmpM));
   
   fprintf(fid, '%s\t%d\n', source.textdata{i+1}, idx);

end

fclose(fid);

display('log/min');

end

function [Vmean] = dsnt(source, target)

    for i = 1:9
        if source(i) > 1
            source(i) = 1;
        elseif source(i) < -1
            source(i) = -1;
        end
    end

    dists = sum((target - repmat(source, [length(target), 1])).^2, 2); 
    Vmean = min(dists);
end




