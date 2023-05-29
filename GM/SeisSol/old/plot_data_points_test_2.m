clear;close all;clc
% Define the file path and file name prefix
file_path = 'data_noVs30_noBuilding';
file_prefix = 'receiver-';

% Define the range of receiver numbers to import
start_num = 1324;
end_num = 1324;

% Loop over the receiver numbers
for i = start_num:end_num
    %% Data import
    file_name = strcat(file_prefix, sprintf('%05d', i), '.dat');
    % Import the data from the file, skipping the first 5 lines
    file_data = importdata(fullfile(file_path, file_name), ' ', 5);
    % Extract the time axis and v(t) data
    time = file_data.data(:, 1);
    v1 = file_data.data(:, end-2);
    v2 = file_data.data(:, end-1);
    v3 = file_data.data(:, end);

    % Find indices corresponding to time range of 1 to 5
    start_idx = find(time >= 1, 1);
    end_idx = find(time <= 5, 1, 'last');

    % Extract the time signal from 1 to 5 and shift it back to 0 to 4
    time = time(start_idx:end_idx) - 1;
    v1 = v1(start_idx:end_idx);
    v2 = v2(start_idx:end_idx);
    v3 = v3(start_idx:end_idx);
    v1_ini=v1;
    v2_ini=v2;
    v3_ini=v3;
    t_ini=time;
    x_label='Time';
    y_label='v(t)';
    fun_plot(time,v1,v2,v3,x_label,y_label,file_name)
    %% Low pass filter
    % Define the filter specifications
    fs = 1 / (time(3) - time(2)); % Sampling frequency in Hz
    f_nqst = fs/2;
    fc = 100; % Cutoff frequency in Hz
    % Compute the normalized cutoff frequency
    fc_norm = fc/fs;

    % Design the filter
    [b, a] = butter(4, fc_norm, 'low');


    % Apply the filter to the data
    v1_filtered = filtfilt(b, a, v1);
    v2_filtered = filtfilt(b, a, v2);
    v3_filtered = filtfilt(b, a, v3);

    x_label='Time';
    y_label='v(t)';
    fun_plot(time,v1_filtered,v2_filtered,v3_filtered,x_label,y_label,file_name)
    u1 = cumtrapz(time, v1_filtered, 1);
    u2 = cumtrapz(time, v2_filtered, 1);
    u3 = cumtrapz(time, v3_filtered, 1);

    %     x_label='Time';
    %     y_label='u(t)';
    %     fun_plot(time,u1,u2,u3,x_label,y_label,file_name)
    t_filt=time;
    v1=v1_filtered;
    v2=v2_filtered;
    v3=v3_filtered;

    %% Resampling
    Fs=200;
    dt=1/Fs;
    t_new=(time(1):dt:time(end)).';
    v1=interp1(time,v1,t_new);
    v2=interp1(time,v2,t_new);
    v3=interp1(time,v3,t_new);
    time=t_new;
    %% FFT
    % Compute the Fourier transform of each v(t) signal
    Fs = 1 / (time(3) - time(2)); % Sampling frequency
    nfft = 2^nextpow2(length(time)); % Number of points for FFT
    freq = Fs / 2 * linspace(0, 1, nfft/2+1); % Frequency axis for FFT
    timestep = 1 / Fs; % Time step
    v1_fft = fft(v1, nfft)* (1/Fs);
    v2_fft = fft(v2, nfft)* (1/Fs);
    v3_fft = fft(v3, nfft)* (1/Fs);
    %% Single sided spectrum
    v1_fft_ss=2*(v1_fft(1:nfft/2+1));
    v2_fft_ss=2*(v2_fft(1:nfft/2+1));
    v3_fft_ss=2*(v3_fft(1:nfft/2+1));
    x_label='f(Hz)';
    y_label='u(f)';
    file_name = strcat(file_prefix, sprintf('%05d', i), '.dat');
    fun_plot(freq,abs(v1_fft_ss),abs(v2_fft_ss),abs(v3_fft_ss),x_label,y_label,file_name)
    %%
    file_name_2 = strcat('fref_vel_f_mat.txt');

    fileID = fopen(fullfile(file_path,file_name_2), 'w');
    % Write the column headers to the file
    fprintf(fileID, 'f(Hz)\tRe\tIm\tAbs\n');

    % Write the data to the file
    for i_f = 1:length(freq)
        fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', freq(i_f), real(v1_fft_ss(i_f)), imag(v2_fft_ss(i_f)), abs(v3_fft_ss(i_f)));
    end
    % Close the file
    fclose(fileID);
    %% Apply zero-padding to the FFT output (single sided spectrum)
    v1_fft_pad=[v1_fft_ss; conj(flipud(v1_fft_ss(1:end-1)))];
    %     v1_fft_pad(end)=[];
    v2_fft_pad=[v2_fft_ss; conj(flipud(v2_fft_ss(1:end-1)))];
    %     v2_fft_pad(end)=[];
    v3_fft_pad=[v3_fft_ss; conj(flipud(v3_fft_ss(1:end-1)))];
    %     v3_fft_pad(end)=[];
    %% Compute the IFFT
    v1_ifft = (0.5)*ifft(v1_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v1_ifft = v1_ifft(1:length(v1));
    v2_ifft = (0.5)*ifft(v2_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v2_ifft = v2_ifft(1:length(v2));
    v3_ifft = (0.5)*ifft(v3_fft_pad*Fs, nfft, 1, 'symmetric'); % Inverse transform with descaling
    v3_ifft = v3_ifft(1:length(v3));

    % Compute the time step and time array
    dt = 1 / Fs;
    t = (0:length(v1_ifft)-1) * dt;
    %% Plotting
    plot_tiles
end