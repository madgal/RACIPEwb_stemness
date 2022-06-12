% min: 3.0476

function ev = entropy_cal_wt(casenum)

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
            tmp(1, i) = 1;
        end
        
    end
    
%     display(sum(tmp));
    
    probs = tmp./sum(tmp);  % with line 19 setted to 1
    
    ev = -1*sum((probs).*log2(probs));

%     ev = tmp./sum(tmp); % with line 19 setted to 0
        
end