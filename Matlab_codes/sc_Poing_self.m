clear;clc;close all
%% Importing Data
evnt='Po2016'

stn_vect={'POI01' 'POI02' 'POI03'};
date='2016_12_20';
time='03_30_51';
r_vect=[1.6 5.6 3.3];

n_stns=length(stn_vect);

bf_nm_u = 'fftd_%d_%s_%s_%s';
bf_nm_v = 'fftv_%d_%s_%s_%s';
cols = {'Freq', 'Re','Im','Amp'};
s_dir = [1 2 3];
n_snr = numel(s_dir);
%%
FACTR_freq = fns_imprtdata.get_FACTR();
FACTR_time = 1;

%%%%
vlbl_vect={'$v_x$,~mm/s' '$v_y$,~mm/s' '$v_z$,~mm/s'};
ulbl_vect={'$u_x$,~mm' '$u_y$,~mm' '$u_z$,~mm'};
%%%%
% figure
for i=1:n_stns
    stn=stn_vect{i}
    ff_fldr = fullfile('GM','GM_POI2016',stn);

    [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    [f_x_max,f_y_max,f_z_max]=fns_plot.find_dominant_freq(ff_Vamp_mat,f_inpt_V{1})
    figure
    % x_lim=[0 80];
    % y_lim=[0 1.4e-4];
    x_lim=[];
    y_lim=[];
    r=r_vect(i);
    fns_plot.plt_ff_svrlstns(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        stn,r,date, time,vlbl_vect,{'Frequency,~Hz'},'cut',evnt,x_lim,y_lim,FACTR_freq)
    if stn=="POI01"
        v_x_Po1=ff_Vamp_mat{1};
        v_y_Po1=ff_Vamp_mat{2};
        v_z_Po1=ff_Vamp_mat{3};
    end
end
% Vinpt_PPV=(v_x_Po1.^2+v_y_Po1.^2+v_z_Po1.^2).^0.5;
% max_Vinpt_PPV=max(Vinpt_PPV)
% [maxrow_indices, ~] = find(bsxfun(@eq, Vinpt_PPV, max_Vinpt_PPV));
% f=f_inpt_V{1};
% f_idx=f(maxrow_indices)
% figure
% plot(f,Vinpt_PPV)
%
% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     ff_fldr = fullfile('GM','GM_POI2016',stn);
%
%     [f_inpt_U,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
%         fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols);
%     % x_lim=[0 80];
%     % y_lim=[0 4e-4];
%     x_lim=[];
%     y_lim=[];
% r=r_vect(i);
%     fns_plot.plt_ff_svrlstns(f_inpt_U, ff_Uamp_mat,bf_nm_u,123,...
%         stn,r,date, time,ulbl_vect,{'f (Hz)'},'cut',evnt,x_lim,y_lim,FACTR_freq)
% end
%%
bf_nm_ut = 'd_%d_%s_%s_%s';
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};
%

for i=1:n_stns
    stn=stn_vect{i}
    ff_fldr = fullfile('GM','GM_POI2016',stn);

    [t_in,ff_Vt]=...
        fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols_t);
    [value_x,value_y,value_z]=fns_plot.find_v_max(ff_Vt)
    % x_lim=[0 80];
    % y_lim=[-0.05 0.05];
    figure
    x_lim=[];
    y_lim=[];
    r=r_vect(i);
    % fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
    %     stn,r,date, time,vlbl_vect,{'t (s)'},'cut',evnt,x_lim,y_lim,FACTR_time)
end

% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     ff_fldr = fullfile('GM','GM_POI2016',stn);
%
%     [t_in,ff_Ut]=...
%         fns_imprtdata.get_ff_tim(bf_nm_ut,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols_t);
%     % x_lim=[0 80];
%     % y_lim=[0 4e-4];
%     x_lim=[];
%     y_lim=[];
% r=r_vect(i);
%     fns_plot.plt_ff_svrlstns(t_in, ff_Ut,bf_nm_ut,123,...
%         stn,r,date, time,ulbl_vect,{'t (s)'},'cut',evnt,x_lim,y_lim,FACTR_time)
% end