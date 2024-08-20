clear;clc;close all
%% Importing Data
% stn_vect={'INSH','TMO54','INS5'};
stn_vect={'TMO54','INS5'};

% stn_vect={'INSH'};

r_vect=[1.9 3 5.2];
date='2013_10_02';
time_evnt='01_13_26';
n_stns=length(stn_vect);


s_dir = [1 2 3];
n_snr = numel(s_dir);
vlbl_vect={'$v_x$ (m/s)' '$v_y$ (m/s)' '$v_z$ (m/s)'};
ulbl_vect={'$u_x$ (m)' '$u_y$ (m)' '$u_z$ (m)'};
%

%%
bf_nm_ut = 'd_%d_%s_%s_%s';
bf_nm_vt = 'v_%d_%s_%s_%s';
bf_nm_at = 'a_%d_%s_%s_%s';
cols_t = {'tim', 'val'};

figure
for i=1:n_stns
    stn=stn_vect{i};
    fldr_nm = ['GM_', stn, '_', date, '_', time_evnt];
    ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time_evnt],...
        fldr_nm);
    [t_in,ff_Vt]=...
        fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
        stn,date, time_evnt,n_snr,ff_fldr,cols_t);
    % fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
    %     stn,date, time_evnt,vlbl_vect,{'t (s)'},'initial',evnt)
end

for i = 1:n_stns
    stn = stn_vect{i};
    fldr_nm = ['GM_', stn, '_', date, '_', time_evnt];
    ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time_evnt], fldr_nm);
    [t_in, Ut] = fns_imprtdata.get_ff_tim(bf_nm_ut, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
    [~, Vt] = fns_imprtdata.get_ff_tim(bf_nm_vt, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
    [~, At] = fns_imprtdata.get_ff_tim(bf_nm_at, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
    for j=1:n_snr
        t_cut=   t_in{j,:};
        mask = (t_cut >= 0) & (t_cut <= 40);
        t_cut = t_cut(mask) - 0;
        t_in{j,:}=t_cut;
        Ut_cut = Ut{:,j}.';
        Ut_cut=Ut_cut(:,mask);
        Utin(:,j)=Ut_cut.';
        Vt_cut = Vt{:,j}.';
        Vt_cut=Vt_cut(:,mask);
        Vtin(:,j)=Vt_cut.';
        At_cut = At{:,j}.';
        At_cut=At_cut(:,mask);
        Atin(:,j)=At_cut.';
    end
    Fs = 1./(t_in{1}(3)-t_in{1}(2));
    ff_fldrnew=fullfile('GM', 'GM_Insheim_Oct2013', stn);
    [u_fft_ss,freq,u_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Utin);
    fns_data_process.save_data_time('d', stn, date, time_evnt,ff_fldrnew,t_cut,Utin);
    fns_data_process.save_data_freq('d', stn, date, time_evnt,ff_fldrnew,freq,u_fft_ss);
    [v_fft_ss,~,v_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Vtin);
    fns_data_process.save_data_time('v', stn, date, time_evnt,ff_fldrnew,t_cut,Vtin);
    fns_data_process.save_data_freq('v', stn, date, time_evnt,ff_fldrnew,freq,v_fft_ss);
    [a_fft_ss,~,a_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Atin);
    fns_data_process.save_data_time('a', stn, date, time_evnt,ff_fldrnew,t_cut,Atin);
    fns_data_process.save_data_freq('a', stn, date, time_evnt,ff_fldrnew,freq,a_fft_ss);
    % Plotting
    leg_vect={'X','Y','Z'};

    x_l_f='Frequency~(Hz)';
    y_l_f='v(f)';
    x_l_t='Time~(s)';
    y_l_t='v(t)';
    plot_tiles.plot_tiles_set_2(2,3,Vtin,abs(v_fft_ss),t_cut,freq,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
end

