function [out, idx, num] = classify_states_Guodata_digit_0(source)
% Manully
% Cdx2	Gata6	Klf4    Nanog	Oct4	Sox2
% type(1,:)   = [-1 -1  1  1  1  1]; 
% type(2,:)   = [ 1  1 -1 -1 -1 -1]; 
% type(3,:)   = [ 1 -1 -1 -1 -1 -1]; 
% type(4,:)   = [-1  1  1 -1  1  1]; 
% type(5,:)   = [-1  1 -1 -1  1 -1];
% type(6,:)   = [-1  1  1  1  1  1];
% type(7,:)   = [ 1 -1 -1  1 -1 -1];
% type(8,:)   = [ 1 -1  1  1  0  1];
% type(9,:)   = [-1 -1  1  1  1  1]; 
% type(10,:)  = [-1 -1  1  1  1 -1];
% type(11,:)  = [ 1 -1  0  0  0  1]; 
% type(12,:)  = [ 1  1  0 -1  0  1];
% type(13,:)  = [ 1  1 -1 -1 -1  1];
% type(14,:)  = [ 1  1 -1  0 -1 -1];
% type(15,:)  = [-1  1  1  1  1 -1]; 
% type(16,:)  = [-1 -1 -1 -1  0 -1];

typetmp0 = load('type.mat');
type = typetmp0.type(:, [2 3 5 6 7 9]);

% type = typetmp(1,:);
% 
% for i = 2:16
%     cnt = 0;
%     for j = 1:(i-1)
%         if ~isequal(typetmp(i,:), typetmp(j,:))
%             cnt = cnt + 1;
%         end
%     end
%     
%     if cnt == (i-1)
%         type = [type; typetmp(i,:)];
%     end
% end

% Automatically
% Cdx2	Gata6	Nanog	Oct4	Sox2
% V = importdata('~/Dropbox/research/rand/stem/c_code/euler/stem_dymT/WT/zscore/wt/Wt_average_expression.dat');
% V2 = zeros(16,9);
% for i = 1:16
%     for j = 1:9
%         if (V.data(i,j) <= 0)
%             V2(i,j) = -1;
%         else
%             V2(i,j) = 1;
%         end
%     end
% end
% 
% for i = 1:16
%     if V2(i, 7) < V2(i, 8)
%         V2(i,7) = V2(i,8);
%     end
% 
%     if V2(i, 9) < V2(i, 8)
%         V2(i,9) = V2(i,8);
%     end
% end
%         
% typetmp = V2(:, [2, 3, 5, 6, 7, 9]);
% % typetmp(9, :)  = [-1 -1 -1 -1  1  1];
% % typetmp(16,:)  = [-1 -1 -1 -1 -1 -1];
% type = typetmp(1,:);
% 
% for i = 2:16
%     cnt = 0;
%     for j = 1:(i-1)
%         if ~isequal(typetmp(i,:), typetmp(j,:))
%             cnt = cnt + 1;
%         end
%     end
%     
%     if cnt == (i-1)
%         type = [type; typetmp(i,:)];
%     end
% end

% type = type(1:10 ,:); % the state bigger than 5%

display(type);

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
    for j = 1:size(type, 1)
        if isequal(type(j,:), tmp)
            idx(i) = j;
            cnt = 1;
        end
    end
    
    if cnt == 0
        idx(i) = -1;
    end
end

num = size(type, 1);

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



