% gata overexpression

function y2 = dyms(te, num)

meanV = [-2.60795426826831,0.994005361302023,-0.309877169274200,-1.86114775683729,-8.50886949962232,3.13482858057809,-0.207682860023183,1.35545354751524,3.68812332690798];
stdlog= [5.40468372769315,4.28389497337265,7.06476202849428,6.58204259362518,7.59711847918126,3.13861573842060,6.29754529608360,4.94682468998410,3.22170559331197];

fparas = importdata('paras_12.dat');
friv   = importdata('RIV_12.dat');

dt = 0.1;

t = 0:dt:te;
y = zeros(te/dt+1,10);

y(1, 10) = 1;
y(1, 1:9) = friv(num, [3, 1, 2, 5, 6, 4, 7, 8, 9] + 3);
p = num2cell(fparas(num, 4:end));

for i = 1:(length(t) - 1)
    
    ytmp = model(t(i), y(i, :), p{:});
    y(i+1, :) = y(i, :) + ytmp.*dt;

end

y2 = zeros(size(y));
for i = 1:9
    y2(:, i) = (log2(y(:, i)) - meanV(i))./stdlog(i);
end

% figure('Color', [1 1 1]);
% c = subplot(2, 1, 1);
% h = plot(t, (y2(:, 1:9)), 'LineWidth', 2);
% set(h(1),'DisplayName','Gata6', 'Color', [1 0 0]);
% set(h(2),'DisplayName','Gcnf',  'Color', [1 1 1]);
% set(h(3),'DisplayName','Cdx2',  'Color', [1 1 0]);
% set(h(4),'DisplayName','Klf4',  'Color', [1 1 1]);
% set(h(5),'DisplayName','Nanog', 'Color', [1 0 1]);
% set(h(6),'DisplayName','Pbx1',  'Color', [1 1 1]);
% set(h(7),'DisplayName','Oct4',  'Color', [0 0 1]);
% set(h(8),'DisplayName','OS',    'Color', [0 1 1]);
% set(h(9),'DisplayName','Sox2',  'Color', [0 1 0]);
% legend1 = legend(c,'show');
% set(legend1,'FontSize',9);
% 
% subplot(2, 1, 2);
% plot(t, y(:, 10), 'LineWidth', 2);

% Map to PCA
% Edge{1} = linspace(-5, 5, 100);
% Edge{2} = linspace(-3, 3, 100);
% 
% V = importdata('../stem_solution_0_zscore.dat');
% 
% [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(V.data, 'Algorithm', 'svd');
% 
% figure;
% hold on;
% hist3(SCORE(:, 1:2), 'Edges',Edge);
% set(get(gca, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto', 'LineStyle','none');
% axis square;
% xlim([-5 5]);
% ylim([-3 3]);
% view(2);
% 
% G = zeros(size(y2, 1), 2);
% idx = [2 3 1 6 4 5 7 8 9];
% 
% for i = 1:9
%     G(:, 1) = G(:, 1) + COEFF(i, 1).*y2(:, idx(i));
%     G(:, 2) = G(:, 2) + COEFF(i, 2).*y2(:, idx(i));
% end
% 
% plot3(G(:, 1), G(:, 2), ones(size(G, 1), 1).*1000, 'wo-');


end


function out = model(t, y, ga, gb, gc, gd, ge, gf, gg, gh, gI, ka, kb, kc, kd, ke, kf, kg, kh, kI, lamdaaa, lamdaha, lamdaea, lamdaga, lamdaab, lamdacb, lamdacc, lamdaec, lamdagc, lamdaId, lamdaed, lamdagd, lamdaee, lamdafe, lamdahe, lamdade, lamdaae, lamdace, lamdaef, lamdabg, lamdacg, lamdahg, lamdaIh,  lamdagh, lamdahI, a0a, h0a, e0a, g0a, a0b, c0b, c0c, e0c, g0c, I0d, e0d, g0d, e0e, f0e, h0e, d0e, a0e, c0e, e0f, b0g, c0g, h0g, I0h, g0h, h0I, naa, nha, nea, nga, nab, ncb, ncc, nec, ngc, nId, ned, ngd, nee, nfe, nhe, nde, nae, nce, nef, nbg, ncg, nhg, nIh, ngh, nhI)

  if t <= 500
      rt =  0.1;
  elseif t> 1000
      rt = 0;
  else
      rt = -0.1;
  end

  out = [ga*Hillshift(y(1), a0a, naa, lamdaaa)*Hillshift(y(8), h0a, nha, lamdaha)*Hillshift(y(5), e0a, nea, lamdaea)*Hillshift(y(7), g0a, nga, lamdaga)  - ka * y(1), ... %A Gata6
         gb*Hillshift(y(1), a0b, nab, lamdaab)*Hillshift(y(3), c0b, ncb, lamdacb)  - kb * y(2), ... %B Gcnf
  	     gc*Hillshift(y(3), c0c, ncc, lamdacc)*Hillshift(y(5), e0c, nec, lamdaec)*Hillshift(y(7), g0c, ngc, lamdagc) - kc * y(3),  ... %C Cdx2
         gd*Hillshift(y(9), I0d, nId, lamdaId)*Hillshift(y(5), e0d, ned, lamdaed)*Hillshift(y(7), g0d, ngd, lamdagd) - kd * y(4),  ... %D Klf4
         ge*Hillshift(y(5), e0e, nee, lamdaee)*Hillshift(y(6), f0e, nfe, lamdafe)*Hillshift(y(8), h0e, nhe, lamdahe)*Hillshift(y(4), d0e, nde, lamdade)*Hillshift(y(1), a0e, nae, lamdaae)*Hillshift(y(3), c0e, nce, lamdace) - ke * y(5), ... %E Nanog
         gf*Hillshift(y(5), e0f, nef, lamdaef) - kf * y(6), ... %F Pbx1
         gg*Hillshift(y(2), b0g, nbg, lamdabg)*Hillshift(y(3), c0g, ncg, lamdacg)*Hillshift(y(8), h0g, nhg, lamdahg) - kg * y(7), ... %G Oct4
         gh*Hillshift(y(9), I0h, nIh, lamdaIh)*Hillshift(y(7), g0h, ngh, lamdagh) - kh * y(8), ... %H Oct4-Sox2
         (gI/y(10))*Hillshift(y(8), h0I, nhI, lamdahI) - kI*y(9), ... %I Sox2
         rt];
     
end

function out = Hillshift(x, x0, nx, lamda)
    out = lamda + (1.0 - lamda) .* (1.0./(1.0 + (x./x0).^nx));
end
