f_vect=(0:1:100).';
% K_vect=(sin(1,length(f_vect))).';
% K_vect=1e9*(ones(1,length(f_vect))).';
% K_vect(2,1)=1e4;

K_vect=triangularPulse(10,30,50,f_vect)*(1e3-1.32E+10)+1.32E+10;
figure
plot(f_vect,K_vect)

data = [f_vect K_vect];

fid = fopen('K_omega_1_2.txt', 'w');
fprintf(fid, '%14.7e %14.7e\n', data');
% sprintf('%u %u \n', [A B]')
fclose(fid);
type 'K_omega_1_2.txt'
