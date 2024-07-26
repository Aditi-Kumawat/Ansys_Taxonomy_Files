%% Initialization
clear;clc;close all

rf_fldr = fns_Inpt_BldPara.selectRfFldr();
disp(['selectedRfFldr: ', rf_fldr])
%----------------------------------------------------%
bld_soil_fndn = fns_Inpt_BldPara.select_bldsoilpara();
disp(['bld_soil_fndnPara: ', bld_soil_fndn])
%----------------------------------------------------%
name_evnt= fns_EvntData.select_event_stn();
disp(['selected_event: ', name_evnt])

[l_vect,b_vect,h,wall_config,dampg_vect,bld_cases]=fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);

[n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);

[evnt,stn,~,date,time,nzero,ff_fldr,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=...
    fns_EvntData.get_event_stn(name_evnt);
%% Importing Data
[f_inpt,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
    fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
    stn,date, time,n_snr,ff_fldr,cols);
[f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
    fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
    stn,date, time,n_snr,ff_fldr,cols);
fns_plot.plt_ff(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
    stn,date, time,'Velocity~(m/s)','initial')
%% Importing Transfer Function
n_c = length(cmpt); % three components
v_ref=5e-8;
for i_str = 0:n_str
    [f_vect, TF_amp_mat, TF_cpmlx_mat] = fns_Wall_and_DR.get_TF(rf_fldr,...
        n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f,...
        wall_config, dampg_vect, i_str, n_c);

    for i_c = 1:n_c
        %% Calculating Velocity
        TFcpmlx_intrp{i_c}=interp1(f_vect,TF_cpmlx_mat{i_c},...
            f_inpt{i_c},'linear','extrap');
        Ucmplx_mat{i_c}=ff_Ucmplx_mat{i_c}.*TFcpmlx_intrp{i_c};

        Vabs_mat{i_c}=abs(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
        Vcmplx_mat{i_c}=(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
        Vdb_mat{i_c}=20*log10(Vabs_mat{i_c}(nzero:end,:)./v_ref);
    end
    %%
    Vdb_zCell{i_str+1} = Vdb_mat{3};
    Vdb_xCell{i_str+1} = Vdb_mat{1};
    Vdb_yCell{i_str+1} = Vdb_mat{2};
    f_dBz=f_inpt_V{3}(nzero:end);
    f_dBy=f_inpt_V{2}(nzero:end);
    f_dBx=f_inpt_V{1}(nzero:end);
end

%%
ylblvectz = {'$$v_{z}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};
ylblvecty = {'$v_{y}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};
ylblvectx = {'$v_{x}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};
ylblvect = {'$v_{h}$,~dB'; '[ref: $5 \times 10^{-8}$ m/s]'};

dfz=f_dBz(3)-f_dBz(2);
dfy=f_dBy(3)-f_dBy(2);
dfx=f_dBx(3)-f_dBx(2);
num_lb=bld_cases;
for i_flur = 1:n_str+1
% for i_flur = 1
    %     [f,v,v_db]=fns_unitgeomdb.DIN4150_3_lims(v_ref,i_flur);
    Vdb_Zmat = Vdb_zCell{i_flur};
    Vdb_Xmat = Vdb_xCell{i_flur};
    Vdb_Ymat = Vdb_yCell{i_flur};
    min_l = min([length(Vdb_Xmat), length(Vdb_Ymat), length(Vdb_Zmat)]);
    Vdb_mat=(Vdb_Zmat(1:min_l,:).^2+Vdb_Xmat(1:min_l,:).^2+Vdb_Ymat(1:min_l,:).^2).^0.5;
    Vdb_horz_mat=(Vdb_Xmat(1:min_l,:).^2+Vdb_Ymat(1:min_l,:).^2).^0.5;
    for j=1:num_lb
        % Vdb_vect=Vdb_mat(:,j);
        % [Fun_rms_vect, f_cenVect] = fns_Octve.get_octBand(...
        %     Vdb_vect, f_dBz, dfz);
        Vdb_horz_vect=Vdb_horz_mat(:,j);
        [Fun_rms_vect_horz, f_cenVect_horz] = fns_Octve.get_octBand(...
            Vdb_horz_vect, f_dBx, dfx);
        V_rms_horz_mat(:, j)=Fun_rms_vect_horz;
        Vdb_Zvect=Vdb_Zmat(:,j);
        [Fun_rms_vectz, f_cenVectz] = fns_Octve.get_octBand(...
            Vdb_Zvect, f_dBz, dfz);
        V_rms_Zmat(:, j)=Fun_rms_vectz;
        Vdb_Yvect=Vdb_Ymat(:,j);
        [Fun_rms_vecty, f_cenVecty] = fns_Octve.get_octBand(...
            Vdb_Yvect, f_dBy, dfy);
        V_rms_Ymat(:, j)=Fun_rms_vecty;
        Vdb_Xvect=Vdb_Xmat(:,j);
        [Fun_rms_vectx, f_cenVectx] = fns_Octve.get_octBand(...
            Vdb_Xvect, f_dBx, dfx);
        V_rms_Xmat(:, j)=Fun_rms_vectx;
    end
    % fns_unitgeomdb.plt_Vrms(f_cenVect,V_rms_Zmat,f_iso,i_flur,V_s,ylblvect{i_flur})
    y_lim=[0 85];
    fns_unitgeomdb.plt_Vrms_stats(f_cenVectz,V_rms_Zmat,i_flur,V_s,ylblvectz,'Z',stn,n_str,y_lim,rf_fldr);
    fns_unitgeomdb.plt_Vrms_stats(f_cenVecty,V_rms_Ymat,i_flur,V_s,ylblvecty,'Y',stn,n_str,y_lim,rf_fldr);
    % fns_unitgeomdb.plt_Vrms_stats(f_cenVect_horz,V_rms_horz_mat,i_flur,V_s,ylblvect,'horz',stn,n_str,y_lim,rf_fldr);
    fns_unitgeomdb.plt_Vrms_stats(f_cenVectx,V_rms_Xmat,i_flur,V_s,ylblvectx,'X',stn,n_str,y_lim,rf_fldr);
    %% plot velocity results in db for individual direction for linear f scale
    %     Vzdb_mean = mean(Vdb_Zmat, 2);
    %     Vzdb_std = std(Vdb_Zmat, 0, 2);
    %     vz_mean = Vzdb_mean;
    %     vz_SDlow = Vzdb_mean - Vzdb_std;
    %     vz_SDup = Vzdb_mean + Vzdb_std;
    %     %     fns_unitgeomdb.plt_Vdb_MeanSD(f_inptz,vz_mean,vz_SDlow,vz_SDup,...
    %     %         0,100,-20,120,i_flur,V_s,rf_fldr,'Z','$v_z$ (dB; ref: 1nm/s)')
    %     Vdb_Xmat = Vdb_xCell{i_flur};
    %     Vxdb_mean = mean(Vdb_Xmat, 2);
    %     Vxdb_std = std(Vdb_Xmat, 0, 2);
    %     vx_mean = Vxdb_mean;
    %     vx_SDlow = Vxdb_mean - Vxdb_std;
    %     vx_SDup = Vxdb_mean + Vxdb_std;
    %     fns_unitgeomdb.plt_Vdb_MeanSD(f_inptx,vx_mean,vx_SDlow,vx_SDup,...
    %         0,100,-20,120,i_flur,V_s,rf_fldr,'X','$v_x$ (dB; ref: 1nm/s)')
end

