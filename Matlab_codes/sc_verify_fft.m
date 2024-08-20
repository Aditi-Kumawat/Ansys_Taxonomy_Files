clear;close all;clc
% ------------------Verification for FFT/IFFT--------------------
dt=0.005;
t_in=0:dt:15;
t_in=t_in.';
Fs=1./dt; % sampling rate
n_fft=length(t_in); %number of samples

%frequnecy array
df=(Fs/length(t_in));
freq = 0:df:Fs;

% Example function (you should replace this with your actual data)
random_fn = sin(2 * pi * 4 * t_in) + 0.5 * sin(2 * pi * 10 * t_in);
% Calculate the Fourier Transform of velocity
random_fn_fft = fft(random_fn)/Fs;  %dividing by sampling rate
% single-sided fft and frequency
random_fn_fft_ss = random_fn_fft(1:n_fft/2+1);
freq_ss = freq(1:n_fft/2+1);

figure
plot(freq_ss,abs(random_fn_fft_ss));
legend('Fourier Transform');  
% ----------------------IFFT---------------------------------------
random_fn_fft_pad = [random_fn_fft_ss; conj(flipud(random_fn_fft_ss(2:end-1,:)))];
random_fn_ifft  = ifft(random_fn_fft_pad*Fs, n_fft, 1, 'symmetric');  %multiplying by sampling rate
% random_fn_ifft = random_fn_ifft(1:length(t_in),:);
figure
plot(t_in,(random_fn));
hold on
plot(t_in,(random_fn_ifft),'-.r');
legend('Original Function', 'Inverse Fourier Transform');  