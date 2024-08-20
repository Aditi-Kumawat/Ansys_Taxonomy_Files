clear;clc;close all
%----------------------------------------------------%
data_set= fns_EvntData.select_event_stn();
disp(['selected_dataset: ', data_set])
[event_vect,nzero,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=...
    fns_EvntData.getEventList(data_set);
%%
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};
%----------------------------------------------------%
vlbl_vect={'$v_x$ (m/s)' '$v_y$ (m/s)' '$v_z$ (m/s)'};
ulbl_vect={'$u_x$ (m)' '$u_y$ (m)' '$u_z$ (m)'};
%----------------------------------------------------%
%% Start a parallel pool
if isempty(gcp('nocreate'))
    parpool;
end
parfor ie = 1:length(event_vect)
    evnt=event_vect{ie};
    disp(['evnt: ', evnt])
    [stn_vect,date_evnt,time_evnt,r_vect]=fns_data_process.get_event_fordataprocess(data_set,evnt);
    n_stns=length(stn_vect);
    % COMMENT/UNCOMMENT "figure" below or inside the for loop depending upon
    % if you want all stns in a single figure file or each stn in a different figure
    % figure
    for i=1:n_stns
        stn=stn_vect{i};
        [ff_fldr]=fns_EvntData.get_ff_fldr1(data_set,stn,date_evnt,time_evnt);
        [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
            fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
            stn,date_evnt, time_evnt,n_snr,ff_fldr,cols);
        fns_CloudAnalysis.save_freq_forVmax(ff_Vamp_mat,f_inpt_V,ff_fldr)
        r=r_vect(i);
        % x_lim=[0 80];
        % y_lim=[0 4e-4];
        x_lim=[];
        % y_lim=[0 max_val];
        y_lim=[];
        % COMMENT/UNCOMMENT "figure" below or above for loop depending upon
        % if you want all stns in a single figure file or each stn in a different figure
        % figure
        % fns_plot.plt_ff_svrlstns(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        %     stn,r,date_evnt, time_evnt,vlbl_vect,{'f (Hz)'},'initial',evnt,x_lim,y_lim,FACTR_freq)
    end

    % COMMENT/UNCOMMENT "figure" below or inside the for loop depending upon
    % if you want all stns in a single figure file or each stn in a different figure
    % figure
    for i=1:n_stns
        stn=stn_vect{i};
        [ff_fldr]=fns_EvntData.get_ff_fldr1(data_set,stn,date_evnt,time_evnt);
        [t_in,ff_Vt]=...
            fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
            stn,date_evnt, time_evnt,n_snr,ff_fldr,cols_t);
        fns_CloudAnalysis.save_PGV(ff_Vt,ff_fldr)
        r=r_vect(i);
        % COMMENT/UNCOMMENT "figure" below or above for loop depending upon
        % if you want all stns in a single figure file or each stn in a different figure
        % figure
        % x_lim=[0 40];
        % y_lim=[0 4e-4];
        x_lim=[];
        y_lim=[];
        % fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
        %     stn,r,date_evnt, time_evnt,vlbl_vect,{'t (s)'},'initial',evnt,x_lim,y_lim,FACTR_time_evnt)
    end
end