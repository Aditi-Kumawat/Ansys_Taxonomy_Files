clear;clc;close all
%% Importing Data
% evnt={'Po2016', 'Po2017'};
evnt='Po2016'
stn_vect={'POI01'};

% stn_vect={'POI01', 'POI02', 'POI03'};
date='2016_12_20';
time='03_30_51';


n_stns=length(stn_vect);

bf_nm_u = 'fftd_%d_%s_%s_%s';
bf_nm_v = 'fftv_%d_%s_%s_%s';
cols = {'Freq', 'Re','Im','Amp'};
s_dir = [1 2 3];
n_snr = numel(s_dir);
%% 
FACTR_freq = fns_imprtdata.get_FACTR();
FACTR_time = 1;
%%%%%
vlbl_vect={'$v_x$ (m/s)' '$v_y$ (m/s)' '$v_z$ (m/s)'};
ulbl_vect={'$u_x$ (m)' '$u_y$ (m)' '$u_z$ (m)'};

figure
for i=1:n_stns
    stn=stn_vect{i};
    fldr_nm = ['GM_', stn, '_', date, '_', time];
    ff_fldr = fullfile('GM', 'GM_Poing_Sync', ['GM_', date, '_', time],...
        fldr_nm);
    [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    % x_lim=[0 80];
    % y_lim=[0 4e-4];
    x_lim=[];
    y_lim=[];
    fns_plot.plt_ff_svrlstns(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        stn,date, time,vlbl_vect,{'f (Hz)'},'initial',evnt,x_lim,y_lim,FACTR_freq)
end

% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     fldr_nm = ['GM_', stn, '_', date, '_', time];
%     ff_fldr = fullfile('GM', 'GM_Poing_Sync', ['GM_', date, '_', time],...
%         fldr_nm);
%     [f_inpt_U,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
%         fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols);
%     % x_lim=[0 80];
%     % y_lim=[0 4e-4];
%     x_lim=[];
%     y_lim=[];
%     fns_plot.plt_ff_svrlstns(f_inpt_U, ff_Uamp_mat,bf_nm_u,123,...
%         stn,date, time,ulbl_vect,{'f (Hz)'},'initial',evnt,x_lim,y_lim,FACTR_freq)
% end
%%
bf_nm_ut = 'd_%d_%s_%s_%s';
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};

figure
for i=1:n_stns
    stn=stn_vect{i};
    fldr_nm = ['GM_', stn, '_', date, '_', time];
    ff_fldr = fullfile('GM', 'GM_Poing_Sync', ['GM_', date, '_', time],...
        fldr_nm);
    [t_in,ff_Vt]=...
        fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols_t);
    % x_lim=[0 80];
    % y_lim=[0 4e-4];
    x_lim=[];
    y_lim=[];
    fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
        stn,date, time,vlbl_vect,{'t (s)'},'initial',evnt,x_lim,y_lim,FACTR_time)
end

% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     fldr_nm = ['GM_', stn, '_', date, '_', time];
%     ff_fldr = fullfile('GM', 'GM_Poing_Sync', ['GM_', date, '_', time],...
%         fldr_nm);
%     [t_in,ff_Ut]=...
%         fns_imprtdata.get_ff_tim(bf_nm_ut,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols_t);
%     % x_lim=[0 80];
%     % y_lim=[0 4e-4];
%     x_lim=[];
%     y_lim=[];
%     fns_plot.plt_ff_svrlstns(t_in, ff_Ut,bf_nm_ut,123,...
%         stn,date, time,ulbl_vect,{'t (s)'},'initial',evnt,x_lim,y_lim,FACTR_time)
% end