function type = guo_state_assgn

type = zeros(9, 6);

for j = 1:9
   V = importdata(['cluster_' num2str(j) '.dat']);
%    H = figure; 
   for i = 1:6; 
       % plot
%        subplot(3,3,i); 
%        histogram(V(:,i), 'Normalization', 'PDF'); 
%        xlim([-4, 4]);
%        hold on; 
%        plot([0 0], [0 1.5], 'r'); 
%        hold off; 

%        digitalize
%          probs = length(find(V(:,i) > 0))/length(V(:,i));
%          if probs > 0.55
%              type(j, i) =  1; 
%          elseif probs <= 0.45
%              type(j, i) = -1;
%          else
%              type(j, i) = 0;
%          end
         
       % digitalize without 0
        Mv = mean(V(:,i));
        if Mv >= 0
            type(j, i) =  1; 
        else
            type(j, i) = -1;
        end
   end 
%    close(H); 

end


end