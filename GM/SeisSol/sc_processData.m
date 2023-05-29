clear;  clc;
close all
% Define file path and prefix
file_path = 'data_Vs30_Bld';
file_prefix = 'receiver-';
x_val = 4000;
y_val = 4000;
% Define receiver range
start_num = 936;
end_num = 936;

% Loop over receiver numbers
for i = start_num:end_num
    %% Data import and extraction
    file_name = sprintf('%s%05d.dat', file_prefix, i);
    file_data = importdata(fullfile(file_path, file_name), ' ', 5);
    time = file_data.data(:, 1);
    v = file_data.data(:, end-2:end);
    time_idx = time >= 1 & time <= 5;
    time = time(time_idx) - 1;
    v = v(time_idx, :);
    v_ini=v;
    t_ini=time;
    %% Filtering
    fs = 1 / (time(3) - time(2));
    fc = 100;
    [b, a] = butter(4, fc/(fs/2), 'low');
    v_filtered = filtfilt(b, a, v);
    % Compute displacement
    u = cumtrapz(time, v_filtered, 1);

    %% Resampling
    Fs = 200;
    t_resampled = (time(1) : 1/Fs : time(end)).';
    v_resampled = interp1(time, v_filtered, t_resampled);
    u_resampled = interp1(time, u, t_resampled);
    time = t_resampled;
    v = v_resampled;
    u = u_resampled;

    %% FFT
    nfft = 2^nextpow2(length(time));
    freq = Fs / 2 * linspace(0, 1, nfft/2+1);
    v_fft = fft(v, nfft) * (1/Fs);
    v_fft_ss = 2 * v_fft(1:nfft/2+1,:);
    v_fft_pad = [v_fft_ss; conj(flipud(v_fft_ss(2:end-1,:)))];
%     v_fft_ss(1:3,:)=0;
    u_fft = fft(u, nfft) * (1/Fs);
    u_fft_ss = 2 * u_fft(1:nfft/2+1,:);
    f_matrix = repmat(freq', 1, 3);
%     u_fft_ss=v_fft_ss./1i./(2*pi*f_matrix);

    u_fft_pad = [u_fft_ss; conj(flipud(u_fft_ss(2:end-1,:)))];

    %% IFFT and plotting
    v_ifft = (0.5)*ifft(v_fft_pad*Fs, nfft, 1, 'symmetric');
    v_ifft = v_ifft(1:length(v(:,1)),:); % Truncate to same length as v
    u_ifft = (0.5)*ifft(u_fft_pad*Fs, nfft, 1, 'symmetric');
    u_ifft = u_ifft(1:length(v(:,1)),:); % Truncate to same length as v
    t = (0:length(v_ifft)-1) / Fs;
    %% Plotting
    leg_vect={'X','Y','Z'};

    x_l_f='Frequency~(Hz)';
    y_l_f='u(f)';
    x_l_t='Time~(s)';
    y_l_t='u(t)';
    %     plot_tiles.plot_tiles_set_3(3,3,v_ini,v_resampled,v_ifft,t_ini,t_resampled,t,x_l_t,y_l_t,leg_vect)
%     plot_tiles.plot_tiles_set_2(2,3,u_resampled,abs(u_fft_ss),t_resampled,freq,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
%     plot_tiles.plot_tiles_set_2(2,3,u_resampled,u_ifft,t_resampled,t,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
    y_l_f='v(t)';
    y_l_t='v(t)';
%     plot_tiles.plot_tiles_set_2(2,3,v_resampled,abs(v_fft_ss),t_resampled,freq,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
    plot_tiles.plot_tiles_set_2(2,3,v_resampled,v_ifft,t_resampled,t,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
    %% Save data for Ansys
    save_data
end

