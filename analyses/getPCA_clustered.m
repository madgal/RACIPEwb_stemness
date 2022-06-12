file = sprintf('clustPCA_SC_og.dat',i);
m = csvread(file);
levels = m(:,3:11);
muO = mean(levels);
sigO= std(levels);
levels =(levels-muO)./sigO;
[coeffOG,score,latent,tsquared,explainedOG,mu] = pca(levels,'NumComponents',2,'Algorithm','svd');
FF = [m(:,1) m(:,2) score];
FF = FF.';
fileout = 'PCAresults_SC_og.dat';
fid = fopen(fileout,'w');
fprintf(fid,'%6.4f,%6.4f,,\n',explainedOG(1:2));
fprintf(fid,'%2d,%2d,%6.4f,%6.4f\n',FF);
fclose(fid);
fileout = 'PCAresults_N_og.dat';
fid = fopen(fileout,'w');
fprintf(fid,'%6.4f,%6.4f,,\n',explainedOG(1:2));
fprintf(fid,'%2d,%2d,%6.4f,%6.4f\n',FF);
fclose(fid);
numDir = 3

for i=1:numDir
	file = sprintf('clustPCA_SC_sf%d.dat',i);
	m = csvread(file);
	levels = m(:,3:11);
	mu = mean(levels);
	sig= std(levels);
	levels =(levels-mu)./sig;
	[coeff,score,latent,tsquared,explained,mu] = pca(levels,'NumComponents',2,'Algorithm','svd');
	FF = [m(:,1) m(:,2) score];
	FF = FF.';
	fileout = sprintf('PCAresults_SC_sf%d.dat',i);
	fid = fopen(fileout,'w');
	fprintf(fid,'%6.4f,%6.4f,,\n',explained(1:2));
	fprintf(fid,'%2d,%2d,%6.4f,%6.4f\n',FF);
	fclose(fid);

	file = sprintf('clustPCA_N_sf%d.dat',i);
	m = csvread(file);
	levels = m(:,3:11);
	%levels = m(:,3:11).*m(:,3);
	levels =(levels-muO)./sigO;
	score = levels*coeffOG;
	FF = [m(:,1) m(:,2) score];
	FF = FF.';
	fileout = sprintf('PCAresults_N_sf%d.dat',i);
	fid = fopen(fileout,'w');
	fprintf(fid,'%6.4f,%6.4f,,\n',explainedOG(1:2));
	fprintf(fid,'%2d,%2d,%6.4f,%6.4f\n',FF);
	fclose(fid);
end
