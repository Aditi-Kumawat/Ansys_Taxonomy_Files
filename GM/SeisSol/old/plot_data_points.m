clear;close all;clc

% Define the file path and file name prefix
file_path = 'data_noVs30_noBuilding';
file_prefix = 'receiver-';

% Define the range of receiver numbers to import
start_num = 1324;
end_num = 1324;

% Loop over the receiver numbers
for i = start_num:end_num
    % Construct the file name

    file_name = strcat(file_prefix, sprintf('%05d', i), '.dat');

    % Import the data from the file, skipping the first 5 lines
    file_data = importdata(fullfile(file_path, file_name), ' ', 5);

    % Extract the time axis and v(t) data
    time = file_data.data(:, 1);
    v1 = file_data.data(:, end-2);
    v2 = file_data.data(:, end-1);
    v3 = file_data.data(:, end);


    %    file_name_1 = strcat(file_prefix, sprintf('%05d_vt', i), '.txt');
    file_name_1 = strcat('fref_vel_t_mat.txt');

    % Open the file for writing
    fileID = fopen(fullfile(file_path,file_name_1), 'w');

    % Write the column headers to the file
    fprintf(fileID, 't(s)\tv_x(t)\tv_y(t)\tv_z(t)\n');

    % Write the data to the file
    for i_v = 1:length(time)
        fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', time(i_v), v1(i_v), v2(i_v), v3(i_v));
    end
    % Close the file
    fclose(fileID);
    % Compute the Fourier transform of each v(t) signal
    Fs = 1 / (time(2) - time(1)); % Sampling frequency
    nfft = 2^nextpow2(length(time)); % Number of points for FFT
    freq = Fs / 2 * linspace(0, 1, nfft/2+1); % Frequency axis for FFT
    timestep = 1 / Fs; % Time step
    v1_fft = fft(v1, nfft)* (1/Fs);
    v2_fft = fft(v2, nfft)* (1/Fs);
    v3_fft = fft(v3, nfft)* (1/Fs);
    %% Single sided spectrum
    v1_fft_ss=(v1_fft(1:nfft/2+1));
    v2_fft_ss=(v2_fft(1:nfft/2+1));
    v3_fft_ss=(v3_fft(1:nfft/2+1));
    %%
    %     freq_100 = freq(freq <= 100);
    %     v1_fft_ss_100=v1_fft_ss(1:length(freq));
    %     v2_fft_ss_100=v2_fft_ss(1:length(freq));
    %     v3_fft_ss_100=v3_fft_ss(1:length(freq));
    for i_v = 1:3
        if i_v == 1
            v_fft_ss = v1_fft_ss;
        elseif i_v == 2
            v_fft_ss = v2_fft_ss;
        else
            v_fft_ss = v3_fft_ss;
        end

        %         file_name_2 = strcat(file_prefix, sprintf('%05d_vf_%d', i, i_v), '.txt');
        file_name_2 = strcat('fref_vel_f_mat.txt');

        fileID = fopen(fullfile(file_path,file_name_2), 'w');
        % Write the column headers to the file
        fprintf(fileID, 'f(Hz)\tRe\tIm\tAbs\n');

        % Write the data to the file
        for i_f = 1:length(freq)
            fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', freq(i_f), real(v_fft_ss(i_f)), imag(v_fft_ss(i_f)), abs(v_fft_ss(i_f)));
        end
        % Close the file
        fclose(fileID);
    end
    %% Apply zero-padding to the FFT output (single sided spectrum)
    v1_fft_pad=[v1_fft_ss; conj(flipud(v1_fft_ss(1:end-1)))];
    %     v1_fft_pad(end)=[];
    v2_fft_pad=[v2_fft_ss; conj(flipud(v2_fft_ss(1:end-1)))];
    %     v2_fft_pad(end)=[];
    v3_fft_pad=[v3_fft_ss; conj(flipud(v3_fft_ss(1:end-1)))];
    %     v3_fft_pad(end)=[];
    %% Compute the IFFT
    v1_ifft = ifft(v1_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v1_ifft = v1_ifft(1:length(v1));
    v2_ifft = ifft(v2_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v2_ifft = v2_ifft(1:length(v2));
    v3_ifft = ifft(v3_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v3_ifft = v3_ifft(1:length(v3));

    % Compute the time step and time array
    dt = 1 / Fs;
    t = (0:length(v1_ifft)-1) * dt;

    %% Plotting
    % Create a new figure for each receiver file
    figure
    % Plot the v(t) data in three subplots
    subplot(3, 2, 1);
    plot(time, v1);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 1 - %s', file_name));

    subplot(3, 2, 3);
    plot(time, v2);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 2 - %s', file_name));

    subplot(3, 2, 5);
    plot(time, v3);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 3 - %s', file_name));

    % Plot the Fourier transform in three subplots
    subplot(3, 2, 2);
    plot(freq, 2*abs(v1_fft_ss));
    xlabel('Frequency');
    ylabel('Magnitude');
    title(sprintf('Sensor 1 - %s (FFT)', file_name));
    xlim([0,50])
    subplot(3, 2, 4);
    plot(freq, 2*abs(v2_fft_ss));
    xlabel('Frequency');
    ylabel('Magnitude');
    title(sprintf('Sensor 2 - %s (FFT)', file_name));
    xlim([0,50])
    subplot(3, 2, 6);
    plot(freq, 2*abs(v3_fft_ss));
    xlabel('Frequency');
    ylabel('Magnitude');
    title(sprintf('Sensor 3 - %s (FFT)', file_name));
    xlim([0,50])
    % Create a new figure for each receiver file
    figure

    % Plot the v(t) data in three subplots
    subplot(3, 2, 1);
    plot(time, v1);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 1 - %s', file_name));
    subplot(3, 2, 3);
    plot(time, v2);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 2 - %s', file_name));
    subplot(3, 2, 5);
    plot(time, v3);
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 3 - %s', file_name));

    % Plot the inverse Fourier transform of v
    subplot(3, 2, 2);
    plot(t, real(v1_ifft));
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 1 - %s (IFFT)', file_name));
    xlim([0,5])
    subplot(3, 2, 4);
    plot(t, real(v2_ifft));
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 2 - %s (IFFT)', file_name));
    xlim([0,5])
    subplot(3, 2, 6);
    plot(t, real(v3_ifft));
    xlabel('Time');
    ylabel('v(t)');
    title(sprintf('Sensor 3 - %s (IFFT)', file_name));
    xlim([0,5])
end
