V = importdata('stem_solution_0_zscore.dat');
% V = importdata('stem_solution_0_raw.dat');
ID = importdata('stem_solution_0_result_wt_min.dat');

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data, 'Algorithm', 'svd');

mP = zeros(15, 9);

for i = 1:15
    mP(i, :) = mean(SCORE(ID.data == i, :));
end

DP = pdist(mP);
DPM = squareform(DP);

mV = zeros(15, 9);

for i = 1:15
    mV(i, :) = mean(V.data(ID.data == i, :));
end

DV = pdist(mV);
DVM = squareform(DV);

DVM2 = zeros(15, 15);
for i = 1:15
    for j = i+1:15
        DVM2(i, j) = mean((mV(i,:) - mV(j, :)).^2); % MSE l2
        
%         X = mV(i, :);
%         Y = mV(j, :);
%         DVM2(i, j) = 1 - sum((Y - X).^2)/sum((Y - mean(Y)).^2); % RSS

%         DVM2(i, j) = mean((mV(i,:) - mV(j, :))); % MSE l1
    end
end



