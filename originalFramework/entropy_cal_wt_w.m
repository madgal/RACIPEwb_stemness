% wt = 1.2169

function ev = entropy_cal_wt_w(casenum)
    V = importdata(['stem_solution_' num2str(casenum) '_result_wt_min.dat']);
    
    tmp = zeros(1,15);
    
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
    
    p = tmp./sum(tmp);
    
    ev = -1*( ...
            (0.8104 - 0.7449)*eachev(p(2)) + (0.8183 - 0.7449)*eachev(p(3)) + (0.7449 - 0.5633)*eachev(p(2) + p(3)) + ...
            (0.8192 - 0.7878)*eachev(p(14))+ (0.8294 - 0.7878)*eachev(p(5)) + (0.7878 - 0.5633)*eachev(p(14)+ p(5)) + ...
            (0.5633 - 0.4665)*eachev(p(2) + p(3) + p(14) + p(5)) + ...
            (0.8451 - 0.7678)*eachev(p(9)) + (0.8016 - 0.7678)*eachev(p(11))+ (0.7678 - 0.6081)*eachev(p(9) + p(11))+ ...
            (0.8214 - 0.7041)*eachev(p(7)) + (0.8340 - 0.7041)*eachev(p(10))+ (0.7041 - 0.6081)*eachev(p(7) + p(10))+ ...
            (0.6081 - 0.4665)*eachev(p(9) + p(11) + p(7) + p(10)) + ...
            0.4665*eachev(p(2) + p(3) + p(14) + p(5) + p(9) + p(11) + p(7) + p(10)) + ...
            (0.8084 - 0.7445)*eachev(p(1)) + (0.8089 - 0.7445)*eachev(p(12)) + (0.7445 - 0.5856)*eachev(p(1) + p(12))+ ...
            (0.8217 - 0.7682)*eachev(p(6)) + (0.8046 - 0.7682)*eachev(p(4)) + (0.7682 - 0.5856)*eachev(p(4) + p(6)) + ...
            (0.5856 - 0.5243)*eachev(p(1) + p(12) + p(4) + p(6)) + ...
            (0.8132 - 0.7070)*eachev(p(13)) + (0.8062 - 0.7070)*eachev(p(8)) + (0.7070 - 0.5948)*eachev(p(13) + p(8)) + ...
            (0.8189 - 0.5948)*eachev(p(15)) + (0.5948 - 0.5243)*eachev(p(15) + p(13) + p(8)) + ...
            0.5243*eachev(p(1) + p(12) + p(4) + p(6) + p(15) + p(13) + p(8)) ...
            ); 

        
end

function out = eachev(probs)
    out = (probs)*log2(probs);
end
