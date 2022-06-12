function [out, idx] = classify_states_Guodata_digit(type, source, num_type)

idx = zeros(1, length(source.data));

for i = 1:length(source.data)
    tmp = source.data(i, :);
%     for j = 1:6
%        if  tmp(j) > 0
%            tmp(j) =  1;
%        else
%            tmp(j) = -1;
%        end
%     end
    cnt = 0;
    for j = 1:num_type
        if isequal(type(j,:), tmp)
            idx(i) = j;
            cnt = 1;
        end
    end
    
    if cnt == 0
        idx(i) = -1;
    end
end

out = length(source.data) - length(find(idx == -1));

end

function out = cmp_state(x, y)
cnt = 1;
for i = 1:length(x)
    if (x(i) ~= 0 && y(i) ~= 0)
        if x(i) ~= y(i)
            cnt = 0;
        end
    end
end

if cnt == 1
    out = 1;
else
    out = 0;
end

end