function probs = probs_dist(casenum)

    V = importdata(['stem_solution_' num2str(casenum) '_result_wt_min.dat']);
    
    tmp = zeros(1, 15);    
    for i = 1:15
        idx = find(V.data == i);
        
        if ~isempty(idx)
            for j = 1:length(idx)
                name = strsplit(V.textdata{idx(j)}, '_');
                tmp(1, i) = tmp(1, i) + 1/str2double(name{1});
            end
        end
        
        if tmp(1, i) == 0
            tmp(1, i) = 0;
        end

    end
        
    probs = tmp./sum(tmp);
        
end

%probs_plot = probs([14, 5, 4, 6, 1, 12, 8, 13, 10, 15, 7, 11, 9, 3, 2]);