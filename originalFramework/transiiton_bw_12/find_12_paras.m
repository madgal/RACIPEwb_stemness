ID = importdata('../stem_solution_0_result_wt_min.dat');
P0 = importdata('stem_parameters_0.dat');
D  = importdata('stem_solution_0_raw.dat');

fidp = fopen('paras_12.dat', 'w');
fidr = fopen('RIV_12.dat', 'w');

P = P0(P0(:, 2) == 2, :);

for i = 1:(length(ID.data)-1)

    C = strsplit(ID.textdata{i}, '_');
    
    if (str2num(C{1}) == 2 && str2num(C{3}) == 1)
        if ID.data(i)*ID.data(i+1) == 2
           fprintf(fidp, '%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', str2num(C{2}), P(str2num(C{2}), 2), ID.data(i),   P(str2num(C{2}), 3:end)); 
           fprintf(fidp, '%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', str2num(C{2}), P(str2num(C{2}), 2), ID.data(i+1), P(str2num(C{2}), 3:end)); 
           fprintf(fidr, '%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', str2num(C{2}), P(str2num(C{2}), 2), ID.data(i),   D.data(i,   3:end));
           fprintf(fidr, '%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', str2num(C{2}), P(str2num(C{2}), 2), ID.data(i+1), D.data(i+1, 3:end));
        end
    end
    
end

fclose(fidp);
fclose(fidr);