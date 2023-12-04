clear;clc;close all
tic
% List of events
event_vect = {'2009', '2010_1', '2010_2', '2012_1', '2012_2', '2013Jan', '2013Feb', '2013Oct', '2013Nov18', '2013Nov21', '2016_1', '2016_2'};

%%
% FACTR_freq = fns_imprtdata.get_FACTR();
FACTR_time = 1;
%%
s_dir = [1 2 3];
n_snr = numel(s_dir);
vlbl_vect={'$v_x$ (m/s)' '$v_y$ (m/s)' '$v_z$ (m/s)'};
ulbl_vect={'$u_x$ (m)' '$u_y$ (m)' '$u_z$ (m)'};

%%
bf_nm_ut = 'd_%d_%s_%s_%s';
bf_nm_vt = 'v_%d_%s_%s_%s';
bf_nm_at = 'a_%d_%s_%s_%s';
cols_t = {'tim', 'val'};
% Assign a unique vector to each event
parfor ie = 1:length(event_vect)
    evnt=event_vect{ie};
    [stn_vect,date,time_evnt,r_vect]=fns_data_process.get_event_fordataprocess(evnt);
    tmin=0;
    tmax=40;
    n_stns=length(stn_vect);
    %%

    for i = 1:n_stns
        stn = stn_vect{i};
        fldr_evnt=['GM_', date, '_', time_evnt];
        fldr_stn_nm = ['GM_', stn, '_', date, '_', time_evnt];
        ff_fldr = fullfile('GM', 'GM_Insheim_Sync',fldr_evnt , fldr_stn_nm);
        [t_in, Ut] = fns_imprtdata.get_ff_tim(bf_nm_ut, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
        [~, Vt] = fns_imprtdata.get_ff_tim(bf_nm_vt, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
        [~, At] = fns_imprtdata.get_ff_tim(bf_nm_at, s_dir, stn, date, time_evnt, n_snr, ff_fldr, cols_t);
        t_vect=t_in{:,1};
        Nt=length(t_vect);
        Utin=zeros(Nt,n_snr);Vtin=zeros(Nt,n_snr);Atin=zeros(Nt,n_snr);
        for j=1:n_snr
            t_cut=   t_in{j,:};
            mask = (t_cut >= tmin) & (t_cut <= tmax);
            t_cut = t_cut(mask) - tmin;
            Nt1=length(t_cut);
            t_in{j,:}=t_cut;
            Ut_cut = Ut{:,j}.';
            Ut_cut=Ut_cut(:,mask);
            Utin(1:Nt1,j)=Ut_cut.';
            Vt_cut = Vt{:,j}.';
            Vt_cut=Vt_cut(:,mask);
            Vtin(1:Nt1,j)=Vt_cut.';
            At_cut = At{:,j}.';
            At_cut=At_cut(:,mask);
            Atin(1:Nt1,j)=At_cut.';
        end
        Utin(Nt1+1:end,:)=[];
        Vtin(Nt1+1:end,:)=[];
        Atin(Nt1+1:end,:)=[];
        Fs = 1./(t_in{1}(3)-t_in{1}(2));

        cd(fullfile('..', 'GM', 'GM_Insheim_1'));
        if ~exist(fldr_evnt, 'dir')
            mkdir(fldr_evnt);
        end
        ff_fldr_evnt=fullfile('GM', 'GM_Insheim_1', fldr_evnt);
        cd(fullfile('..', '..',ff_fldr_evnt));
        if ~exist(fldr_stn_nm, 'dir')
            mkdir(fldr_stn_nm);
        end
        cd(fullfile('..', '..', '..', 'Matlab_codes'));

        ff_fldrnew=fullfile('GM', 'GM_Insheim_1', fldr_evnt, fldr_stn_nm);
        [u_fft_ss,freq,u_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Utin);
        % fns_data_process.save_data_time('d', stn, date, time_evnt,ff_fldrnew,t_cut,Utin);
        % fns_data_process.save_data_freq('d', stn, date, time_evnt,ff_fldrnew,freq,u_fft_ss);
        [v_fft_ss,~,v_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Vtin);
        % fns_data_process.save_data_time('v', stn, date, time_evnt,ff_fldrnew,t_cut,Vtin);
        % fns_data_process.save_data_freq('v', stn, date, time_evnt,ff_fldrnew,freq,v_fft_ss);
        [a_fft_ss,~,a_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Atin);
        % fns_data_process.save_data_time('a', stn, date, time_evnt,ff_fldrnew,t_cut,Atin);
        % fns_data_process.save_data_freq('a', stn, date, time_evnt,ff_fldrnew,freq,a_fft_ss);
        % Plotting
        leg_vect={'X','Y','Z'};

        x_l_f='Frequency~(Hz)';
        y_l_f='v(f)';
        x_l_t='Time~(s)';
        y_l_t='v(t)';
        % plot_tiles.plot_tiles_set_2(2,3,Vtin,abs(v_fft_ss),t_cut,freq,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
    end
end

processing_time=toc;
disp('total processing time=',processing_time)


