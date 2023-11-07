clear;clc;
%% Importing Data
evnt='Po2016'

stn_vect={'POI01'};
date='2016_12_20';
time='03_30_51';
nzero=6;

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
vlbl_vect={'$v_x$,~m/s' '$v_y$,~m/s' '$v_z$,~m/s'};
ulbl_vect={'$u_x$,~m' '$u_y$,~m' '$u_z$,~m'};
%%%%
v_ref=5e-8;
figure
for i=1:n_stns
    stn=stn_vect{i};
    ff_fldr = fullfile('GM','GM_POI2016',stn);

    [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    % x_lim=[0 80];
    % y_lim=[0 4e-4];
    x_lim=[];
    y_lim=[];

    fns_plot.plt_ff_svrlstns(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        stn,date, time,vlbl_vect,{'Frequency,~Hz'},'cut',evnt,x_lim,y_lim,FACTR_freq)
    if stn=="POI01"
        v_x_Po1=ff_Vamp_mat{1};
        v_y_Po1=ff_Vamp_mat{2};
        v_z_Po1=ff_Vamp_mat{2};
    end
end

f_dBz=f_inpt_V{3}(nzero:end);
dfz=f_dBz(3)-f_dBz(2);
for i_c = 1:3
    Vdb_mat(:,i_c)=20*log10(ff_Vamp_mat{i_c}(nzero:end,:)./v_ref);
    [Fun_rms_vect, f_cenVect] = fns_Octve.get_octBand(...
        Vdb_mat(:,i_c), f_dBz, dfz);
    V_rms_mat(:, i_c)=Fun_rms_vect;
end

figure
subplot(3,1,1);
plot(f_inpt_V{1}(nzero:end,:),Vdb_mat(:,1));
hold on
plot(f_cenVect,V_rms_mat(:,1));
subplot(3,1,2);
plot(f_inpt_V{1}(nzero:end,:),Vdb_mat(:,2));
hold on
plot(f_cenVect,V_rms_mat(:,2));

subplot(3,1,3);
plot(f_inpt_V{1}(nzero:end,:),Vdb_mat(:,3));
hold on
plot(f_cenVect,V_rms_mat(:,3));
y_lim=[0 100];
ylblvectz = {'$v_{z}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};
ylblvectx = {'$v_{x}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};
fns_unitgeomdb.plt_Vrms_stats(f_cenVect,V_rms_mat(:,1),1,450,ylblvectx,'X',stn,0,y_lim)
%%
% bf_nm_ut = 'd_%d_%s_%s_%s';
% bf_nm_vt = 'v_%d_%s_%s_%s';
% cols_t = {'tim', 'val'};
% %
% figure
% for i=1:n_stns
%     stn=stn_vect{i};
%     ff_fldr = fullfile('GM','GM_POI2016',stn);
%
%     [t_in,ff_Vt]=...
%         fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
%         stn,date, time,n_snr,ff_fldr,cols_t);
%     % x_lim=[0 80];
%     % y_lim=[0 4e-4];
%     x_lim=[];
%     y_lim=[];
%     fns_plot.plt_ff_svrlstns(t_in, ff_Vt,bf_nm_vt,123,...
%         stn,date, time,vlbl_vect,{'t (s)'},'cut',evnt,x_lim,y_lim,FACTR_time)
% end