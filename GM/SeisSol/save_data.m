
%    file_name_1 = strcat(file_prefix, sprintf('%05d_vt', i), '.txt');
file_name_1 = strcat(sprintf('vt_rec_%d',i), '.txt');

% Open the file for writing
fileID = fopen(fullfile(file_path,file_name_1), 'w');

% Write the column headers to the file
fprintf(fileID, 't(s)\tv_x(t)\tv_y(t)\tv_z(t)\n');

% Write the data to the file
for i_t = 1:length(time)
    fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', t_resampled(i_t), v_resampled(i_t,1), v_resampled(i_t,2), v_resampled(i_t,3));
end
% Close the file
fclose(fileID);

%%
% file_name_1 = strcat('fref_disp_t_mat.txt');
file_name_1 =strcat(sprintf('ut_rec_%d',i), '.txt');
% Open the file for writing
fileID = fopen(fullfile(file_path,file_name_1), 'w');

% Write the column headers to the file
fprintf(fileID, 't(s)\tu_x(t)\tu_y(t)\tu_z(t)\n');

% Write the data to the file
for i_t = 1:length(time)
    fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', t_resampled(i_t), u_resampled(i_t,1), u_resampled(i_t,2), u_resampled(i_t,3));
end
% Close the file
fclose(fileID);

%%
for i_v=1:3
    file_name_2 = strcat(sprintf('fftv_%d_rec_%d', i_v,i), '.txt');
    fileID = fopen(fullfile(file_path,file_name_2), 'w');

    % Add the lines at the top of the file with the dynamic variables
    fprintf(fileID, 'Spectra of the velocity (only positive half of the frequency range, doubled in amplitude) \n');
    fprintf(fileID, 'Units as in the International System of Units (Hz, m, m/s, m/s^2) \n');
    fprintf(fileID, 'Receiver = %f - x= %f - y= %f  - Channel = %d \n', i, x_val, y_val, i_v);
    fprintf(fileID, '-------------------- \n');

    % Write the column headers to the file
    fprintf(fileID, 'f(Hz)\tRe\tIm\tAbs\n');

    % Write the data to the file
    for i_f = 1:length(freq)
        fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', freq(i_f), real(v_fft_ss(i_f,i_v)), imag(v_fft_ss(i_f,i_v)), abs(v_fft_ss(i_f,i_v)));
    end

    % Close the file
    fclose(fileID);
end


for i_v=1:3
    file_name_2 = strcat(sprintf('fftd_%d_rec_%d', i_v,i), '.txt');
    fileID = fopen(fullfile(file_path,file_name_2), 'w');

    % Add the lines at the top of the file with the dynamic variables
    fprintf(fileID, 'Spectra of the velocity (only positive half of the frequency range, doubled in amplitude) \n');
    fprintf(fileID, 'Units as in the International System of Units (Hz, m, m/s, m/s^2) \n');
    fprintf(fileID, 'Receiver = %f - x= %f - y= %f  - Channel = %d \n', i, x_val, y_val, i_v);
    fprintf(fileID, '-------------------- \n');

    % Write the column headers to the file
    fprintf(fileID, 'f(Hz)\tRe\tIm\tAbs\n');

    % Write the data to the file
    for i_f = 1:length(freq)
        fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', freq(i_f), real(u_fft_ss(i_f,i_v)), imag(u_fft_ss(i_f,i_v)), abs(u_fft_ss(i_f,i_v)));
    end

    % Close the file
    fclose(fileID);
end