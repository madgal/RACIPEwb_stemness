function type = avg_expression

type = zeros(9, 6);

for i = 1:9
    
    V = importdata(['cluster_' num2str(i) '.dat']);
    
% according to mean value 
    Vm = mean(V); 
    idx = find(Vm <= 0);
    type(i, idx) = -1;
    idx2= find(Vm >  0);
    type(i, idx2) = 1;

% according to distribution
%     for j = 1:6
%         
%        probs = length(find(V(:,j) > 0))/length(V(:,j)); 
% 
%        if probs > 0.6
%            type(i, j) = 1;
%        elseif probs <= 0.4
%            type(i, j) = -1;
%        else
%            type(i, j) = 0;
%        end
%            
%     end
end

