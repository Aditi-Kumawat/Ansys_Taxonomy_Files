clear;clc;
t_in=0:0.005:3;
nfft = 2^nextpow2(length(t_in));
Fs = nfft/(t_in(end)-t_in(1));
t_in = transpose(linspace(t_in(1),t_in(end),nfft));
%frequnecy array
freq = Fs / 2 * linspace(0, 1, nfft/2+1);
dfz=freq(3)-freq(2);

%% Verification ----------------------------------------------------------------
% Example velocity function
velocity = sin(2 * pi * 4 * t_in) + 0.5 * sin(2 * pi * 10 * t_in);
% Calculate the Fourier Transform of velocity
velocity_fft = fft(velocity,nfft)*Fs;
velocity_fft_trun = velocity_fft(1:nfft/2+1);
%% IFFT
velocity_fft_pad = [velocity_fft_trun ; conj(flipud(velocity_fft_trun(2:end-1,:)))];
velocity_ifft  = ifft(velocity_fft_pad, nfft, 1, 'symmetric');
velocity_ifft = velocity_ifft(1:length(t_in),:)* (1/Fs);

figure
plot(t_in,velocity);
hold on
plot(t_in,velocity_ifft,'-.');

%%
dt=0.005;
t_in=0:dt:3;
% Example velocity function
vel_g = sin(2 * pi * 4 * t_in) + 0.5 * sin(2 * pi * 10 * t_in);
sample_rate=1/dt;
Nsamples=length(t_in);
df=sample_rate./Nsamples;
freq=(0:df:sample_rate)';

fft_unscaled=fft(vel_g);
% fft scaled
vel_g_bar=fft_unscaled.*(1./sample_rate);
% Back into the time
vel_g_ifft=ifft(vel_g_bar,'symmetric').*sample_rate;

figure
plot(t_in,vel_g);
hold on
plot(t_in,vel_g_ifft,'-.');