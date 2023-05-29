%% Plotting
function fun_plot(data_x,data_1,data_2,data_3,x_label,y_label,file_name)
figure
% Plot the v(t) data in three subplots
subplot(3, 1, 1);
plot(data_x, data_1);
xlabel(x_label);
ylabel(y_label);
title(sprintf('Sensor 1 - %s', file_name));

subplot(3, 1, 2);
plot(data_x, data_2);
xlabel(x_label);
ylabel(y_label);
title(sprintf('Sensor 2 - %s', file_name));

subplot(3, 1, 3);
plot(data_x, data_3);
xlabel(x_label);
ylabel(y_label);
title(sprintf('Sensor 3 - %s', file_name));
end
%     % Plot the Fourier transform in three subplots
%     subplot(3, 2, 2);
%     plot(freq, 2*abs(v1_fft_ss));
%     xlabel('Frequency');
%     ylabel('Magnitude');
%     title(sprintf('Sensor 1 - %s (FFT)', file_name));
%
%     subplot(3, 2, 4);
%     plot(freq, 2*abs(v2_fft_ss));
%     xlabel('Frequency');
%     ylabel('Magnitude');
%     title(sprintf('Sensor 2 - %s (FFT)', file_name));
%
%     subplot(3, 2, 6);
%     plot(freq, 2*abs(v3_fft_ss));
%     xlabel('Frequency');
%     ylabel('Magnitude');
%     title(sprintf('Sensor 3 - %s (FFT)', file_name));
%
%     % Create a new figure for each receiver file
%     figure
%
%     % Plot the v(t) data in three subplots
%     subplot(3, 2, 1);
%     plot(time, v1);
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 1 - %s', file_name));
%     subplot(3, 2, 3);
%     plot(time, v2);
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 2 - %s', file_name));
%     subplot(3, 2, 5);
%     plot(time, v3);
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 3 - %s', file_name));
%
%     % Plot the inverse Fourier transform of v
%     subplot(3, 2, 2);
%     plot(t, real(v1_ifft));
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 1 - %s (IFFT)', file_name));
%     xlim([0,5])
%     subplot(3, 2, 4);
%     plot(t, real(v2_ifft));
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 2 - %s (IFFT)', file_name));
%     xlim([0,5])
%     subplot(3, 2, 6);
%     plot(t, real(v3_ifft));
%     xlabel('Time');
%     ylabel('v(t)');
%     title(sprintf('Sensor 3 - %s (IFFT)', file_name));
%     xlim([0,5])