% 0.3785

cdt = importdata('stem_solution_0_zscore_cdt.dat');
len0 = size(cdt.data, 1);

X = [];
Y = [];
cnt = 0;

len = zeros(1, 15);

for j = 1:15
   V = importdata(['Wt_cluster_' num2str(j) '.dat']);

   tmp = size(V, 1)/len0;
   
   len(j) = size(V, 1);

   if tmp >= 0.005
       cnt = cnt + 1;
       X = [X; V];
       Y = [Y; ones(size(V, 1), 1).*cnt];
   end
end

figure;
[S, h] = silhouette(X, Y);

S2 = zeros(1, 15);
for j = 1:15
    if j == 1
        S2(1) = mean(S(1:len(1)));
    else
        S2(j) = mean(S((sum(len(1:j-1))+1):sum(len(1:j))));
    end
end 
    