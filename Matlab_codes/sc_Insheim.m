clear;clc;close all
%% Importing Data
% evnt='2009'
% evnt='2010_1'
% evnt='2010_2'
% evnt='2012_1'
% evnt='2012_2'
% evnt='2013Jan'
evnt='2013Feb'
% evnt='2013Oct'
% evnt='2013Nov18'
% evnt='2013Nov21'
% evnt='2016_1'
% evnt='2016_2'



rf_fldr = 'input_Data';
if strcmp(evnt, '2009')
    stn_vect={'LDAU'};
    date='2009_10_18';
    time='19_12_12';
elseif strcmp(evnt, '2010_1')
    stn_vect={'LDAU'};
    date='2010_04_07';
    time='09_04_04';
elseif strcmp(evnt, '2010_2')
    stn_vect={'LDAU'};
    date='2010_04_07';
    time='13_46_21';

elseif strcmp(evnt, '2012_1') 
    stn_vect={'INS2','INS3','INS4','INS5','INS6', 'INSH'};
    date='2012_11_12';
    time='11_15_04';
elseif strcmp(evnt, '2012_2') 
    stn_vect={'INS2','INS3','INS4','INS5','INS6', 'INSH'};
    date='2012_11_12';
    time='12_53_02';
elseif strcmp(evnt, '2013Jan') 
    stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7', 'INSH'};
    date='2013_01_26';
    time='19_48_27';
    elseif strcmp(evnt, '2013Oct') %% higest magnitude
    stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
        'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
    date='2013_10_02';
    time='01_13_26';
elseif strcmp(evnt, '2013Feb') %%
    stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7','INSH'};
    date='2013_02_17';
    time='20_07_15';
elseif strcmp(evnt, '2013Nov18') 
    stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
        'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
    date='2013_11_18';
    time='12_54_15';
elseif strcmp(evnt, '2013Nov21') %%
    stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
        'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
    date='2013_11_21';
    time='14_15_22';
elseif strcmp(evnt, '2016_1') 
    stn_vect={'A127A','INS3','INS4B','INS5','INS6B','INS7','INS8',...
        'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55', 'TMO57', 'TMO66'};
    date='2016_02_12';
    time='06_26_04';
elseif strcmp(evnt, '2016_2') 
    stn_vect={'INS1','INS3','INS4B','INS5','INS6B','INS7','INS8',...
        'INSH', 'TMO20', 'TMO54', 'TMO55', 'TMO57', 'TMO66'};
    date='2016_07_14';
    time='17_49_10';
end


n_stns=length(stn_vect);

bf_nm_u = 'fftd_%d_%s_%s_%s';
bf_nm_v = 'fftv_%d_%s_%s_%s';
cols = {'Freq', 'Re','Im','Amp'};
s_dir = [1 2 3];
n_snr = numel(s_dir);
vlbl_vect={'$v_x$ (m/s)' '$v_y$ (m/s)' '$v_z$ (m/s)'};
ulbl_vect={'$u_x$ (m)' '$u_y$ (m)' '$u_z$ (m)'};

figure
for i=1:n_stns
    stn=stn_vect{i};
    fldr_nm = ['GM_', stn, '_', date, '_', time];
    ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time],...
        fldr_nm);
    [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    % x_lim=[0 80];
    % y_lim=[0 4e-4];
    x_lim=[];
    y_lim=[];
    fns_plot.plt_ff_svrlstns(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        stn,date, time,vlbl_vect,{'f (Hz)'},'initial',evnt,x_lim,y_lim)
end

% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     fldr_nm = ['GM_', stn, '_', date, '_', time];
%     ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time],...
%         fldr_nm);
%     [f_inpt_U,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
%         fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols);
%     x_lim=[0 80];
%     y_lim=[0 4e-4];
% x_lim=[];
% y_lim=[];
%     fns_plot.plt_ff_svrlstns(f_inpt_U, ff_Uamp_mat,bf_nm_u,123,...
%         stn,date, time,ulbl_vect,{'f (Hz)'},'initial',evnt,x_lim,y_lim)
% end
% %%
% bf_nm_ut = 'd_%d_%s_%s_%s';
% bf_nm_vt = 'v_%d_%s_%s_%s';
% cols_t = {'tim', 'val'};
%
% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     fldr_nm = ['GM_', stn, '_', date, '_', time];
%     ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time],...
%         fldr_nm);
%     [t_in,ff_Vt]=...
%         fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols_t);
%     x_lim=[0 80];
%     y_lim=[0 4e-4];
% x_lim=[];
% y_lim=[];
%     fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
%         stn,date, time,vlbl_vect,{'t (s)'},'initial',evnt,x_lim,y_lim)
% end
%
% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     fldr_nm = ['GM_', stn, '_', date, '_', time];
%     ff_fldr = fullfile('GM', 'GM_Insheim_Sync', ['GM_', date, '_', time],...
%         fldr_nm);
%     [t_in,ff_Ut]=...
%         fns_imprtdata.get_ff_tim(bf_nm_ut,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols_t);
%     x_lim=[0 80];
%     y_lim=[0 4e-4];
% x_lim=[];
% y_lim=[];
%     fns_plot.plt_ff_svrlstns(t_in, ff_Ut,bf_nm_ut,123,...
%         stn,date, time,ulbl_vect,{'t (s)'},'initial',evnt,x_lim,y_lim)
% end