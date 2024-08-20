clear;clc;close all
tic
% List of events
data_set= fns_EvntData.select_event_stn();
disp(['selected_dataset: ', data_set])
%----------------------------------------------------%
[event_vect,~,~,~,~,~,~,~]=fns_EvntData.getEventList(data_set);
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
for ie = 1:length(event_vect)
    evnt=event_vect{ie};
    [stn_vect,date_evnt,time_evnt,r_vect]=fns_data_process.get_event_fordataprocess(data_set,evnt);
    [tmin,tmax]=fns_data_process.get_tmin_tmax(data_set,evnt);
    n_stns=length(stn_vect);
    %%

    parfor i = 1:n_stns
        stn = stn_vect{i};
        [ff_fldr,fldr_evnt,fldr_stn_nm]=fns_data_process.get_ff_fldr_sync(data_set,stn,date_evnt,time_evnt);

        [t_in, Ut] = fns_imprtdata.get_ff_tim(bf_nm_ut, s_dir, stn, date_evnt, time_evnt, n_snr, ff_fldr, cols_t);
        [~, Vt] = fns_imprtdata.get_ff_tim(bf_nm_vt, s_dir, stn, date_evnt, time_evnt, n_snr, ff_fldr, cols_t);
        [~, At] = fns_imprtdata.get_ff_tim(bf_nm_at, s_dir, stn, date_evnt, time_evnt, n_snr, ff_fldr, cols_t);
        t_vect=t_in{:,1};
        % tmax=t_vect(end);
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

        ff_fldrnew=fns_data_process.create_new_ff_fldr(data_set,fldr_evnt,fldr_stn_nm);

        [u_fft_ss,freq,u_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Utin);
        fns_data_process.save_data_time('d', stn, date_evnt, time_evnt,ff_fldrnew,t_cut,Utin);
        fns_data_process.save_data_freq('d', stn, date_evnt, time_evnt,ff_fldrnew,freq,u_fft_ss);
        [v_fft_ss,~,v_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Vtin);
        fns_data_process.save_data_time('v', stn, date_evnt, time_evnt,ff_fldrnew,t_cut,Vtin);
        fns_data_process.save_data_freq('v', stn, date_evnt, time_evnt,ff_fldrnew,freq,v_fft_ss);
        [a_fft_ss,~,a_ifft]=fns_data_process.fun_fftandifft(t_cut,Fs,Atin);
        fns_data_process.save_data_time('a', stn, date_evnt, time_evnt,ff_fldrnew,t_cut,Atin);
        fns_data_process.save_data_freq('a', stn, date_evnt, time_evnt,ff_fldrnew,freq,a_fft_ss);
        % Plotting
        % leg_vect={'X','Y','Z'};
        % 
        % x_l_f='Frequency~(Hz)';
        % y_l_f='v(f)';
        % x_l_t='Time~(s)';
        % y_l_t='v(t)';
        % plot_tiles.plot_tiles_set_2(2,3,Vtin,abs(v_fft_ss),t_cut,freq,x_l_f,y_l_f,x_l_t,y_l_t,leg_vect)
    end
end



