% average linkage

function classify_states_wt_avg(casenum, num)

name = 'stem_solution';

source = importdata([name '_' num2str(casenum) '_zscore.dat']);

cmpdata = cell(1, num);

for i = 1:num
    cmpdata{i}  = importdata(['Wt_cluster_' num2str(i) '.dat']);
end

fid = fopen([name '_' num2str(casenum) '_result_wt_avg.dat'], 'w');

for i = 1:num
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
   tmpM   = zeros(1,num);
   
   for j = 1:num
       [tmpM(j)] = dsnt(source.data(i,:), cmpdata{j});  
   end
   
   idx = find(tmpM == min(tmpM));
   fprintf(fid, '%s\t%d\n', source.textdata{i+1}, idx);
   
end

fclose(fid);

display('error2_CIE_classification');

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
    Vmean = mean(dists);
end




