%% Initialization
clear;clc;close all;

rf_fldr = fns_Inpt_BldPara.selectRfFldr();
disp(['selectedRfFldr: ', rf_fldr])
%----------------------------------------------------%
bld_soil_fndn = fns_Inpt_BldPara.select_bldsoilpara();
disp(['bld_soil_fndnPara: ', bld_soil_fndn])
%----------------------------------------------------%
name_evnt= fns_EvntData.select_event_stn();
disp(['selected_event: ', name_evnt])

[l_vect,b_vect,h,wall_config,dampg_vect,bldcases]=fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);

[n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);

[evnt,~,stn_vect,date,time,nzero,ff_fldr,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=...
    fns_EvntData.get_event_stn(name_evnt);
%%
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};

n_stns=length(stn_vect);

for i_stn=1:n_stns
    stn=stn_vect{i_stn};
    disp(['stn: ', stn])
    [ff_fldr]=fns_EvntData.get_ff_fldr(name_evnt,stn);
    %% Importing sensor data
    [f_inpt,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
        fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols);
    % fns_plot.plt_ff(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
    %     stn,date, time,'Velocity~(m/s)','initial')

    [t_in,ff_Vt]=...
        fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
        stn,date, time,n_snr,ff_fldr,cols_t);
    t_in=t_in{1};

    n_c = length(cmpt);
    v_ref=5e-8;

    %% Importing Transfer function computed from ANSYS
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
            %% Apply filter transfer to KB(f)
            KB_freq{i_c} = fns_KBvalue.fn_kb_highpass(f_inpt{i_c},Vcmplx_mat{i_c},5.6);
        end
        % Original data
        Vss_zCell{i_str+1} = Vcmplx_mat{3};
        Vss_xCell{i_str+1} = Vcmplx_mat{1};
        Vss_yCell{i_str+1} = Vcmplx_mat{2};

        %Data with filter: KB(f)
        KB_freq_zCell{i_str+1} = KB_freq{3};
        KB_freq_xCell{i_str+1} = KB_freq{1};
        KB_freq_yCell{i_str+1} = KB_freq{2};
    end

    %% IFFT compute v(t)
    Fs = 200;           %sampling rate
    nfft = length(t_in);
    %frequnecy array
    freq = Fs / 2 * linspace(0, 1, nfft/2+1); %length = nfft/2+1
    dfz=freq(3)-freq(2);
    %% ------------------Verification for FFT/IFFT--------------------
    % % Example velocity function (you should replace this with your actual velocity data)
    % random_fn = sin(2 * pi * 4 * t_in) + 0.5 * sin(2 * pi * 10 * t_in);
    % % Calculate the Fourier Transform of velocity
    % random_fn_fft = fft(random_fn)/Fs; %dividing by sampling rate
    % random_fn_fft_trun = random_fn_fft(1:nfft/2+1);
    % figure
    % plot(freq,abs(random_fn_fft_trun));
    % % ----------------------IFFT---------------------------------------
    % random_fn_fft_pad = [random_fn_fft_trun ; conj(flipud(random_fn_fft_trun(2:end-1,:)))];
    % random_fn_ifft  = ifft(random_fn_fft_pad*Fs, nfft, 1, 'symmetric'); %multiplying by sampling rate
    % random_fn_ifft = random_fn_ifft(1:length(t_in),:);
    % figure
    % plot(t_in,(random_fn));
    % hold on
    % plot(t_in,(random_fn_ifft),'-.r');
    %% ----------------------------------------------------------------
    num_lb=bldcases;
    % for i_flur = 1:n_str+1 !!!!!!!! check
    for i_flur = 4
        %% IFFT original data
        Vss_Zmat = Vss_zCell{i_flur};  %length = 2049
        Vss_Xmat = Vss_xCell{i_flur};
        Vss_Ymat = Vss_yCell{i_flur};

        Vzss_fft_pad = [Vss_Zmat; conj(flipud(Vss_Zmat(2:end-1,:)))];
        Vxss_fft_pad = [Vss_Xmat; conj(flipud(Vss_Xmat(2:end-1,:)))];
        Vyss_fft_pad = [Vss_Ymat; conj(flipud(Vss_Ymat(2:end-1,:)))];

        Vz_ifft = ifft(Vzss_fft_pad*Fs, nfft, 1, 'symmetric');
        Vx_ifft = ifft(Vxss_fft_pad*Fs, nfft, 1, 'symmetric');
        Vy_ifft = ifft(Vyss_fft_pad*Fs, nfft, 1, 'symmetric');
        Vz_ifft = Vz_ifft(1:length(t_in),:);
        Vx_ifft = Vx_ifft(1:length(t_in),:);
        Vy_ifft = Vy_ifft(1:length(t_in),:);

        % Store the result in Cell : {floor num}
        Vz_ifft_cell{i_flur} = Vz_ifft;
        Vx_ifft_cell{i_flur} = Vx_ifft;
        Vy_ifft_cell{i_flur} = Vy_ifft;


        %% IFFT data with filter
        KBss_Zmat = KB_freq_zCell{i_flur};
        KBss_Xmat = KB_freq_xCell{i_flur};
        KBss_Ymat = KB_freq_yCell{i_flur};
        KBz_fft_pad = [KBss_Zmat; conj(flipud(KBss_Zmat(2:end-1,:)))];
        KBx_fft_pad = [KBss_Xmat; conj(flipud(KBss_Xmat(2:end-1,:)))];
        KBy_fft_pad = [KBss_Ymat; conj(flipud(KBss_Ymat(2:end-1,:)))];

        KBz_ifft = ifft(KBz_fft_pad*Fs, nfft, 1, 'symmetric');
        KBx_ifft = ifft(KBx_fft_pad*Fs, nfft, 1, 'symmetric');
        KBy_ifft = ifft(KBy_fft_pad*Fs, nfft, 1, 'symmetric');
        KBz_ifft = KBz_ifft(1:length(t_in),:);
        KBx_ifft = KBx_ifft(1:length(t_in),:);
        KBy_ifft = KBy_ifft(1:length(t_in),:);

        % Store the result in Cell : {floor num}
        KBz_ifft_cell{i_flur} = KBz_ifft;
        KBx_ifft_cell{i_flur} = KBx_ifft;
        KBy_ifft_cell{i_flur} = KBy_ifft;

        t = (0:length(Vz_ifft)-1) / Fs;

        max_Vxmat(i_stn,:,i_flur)=max(Vx_ifft);
        max_Vymat(i_stn,:,i_flur)=max(Vy_ifft);
        max_Vzmat(i_stn,:,i_flur)=max(Vz_ifft);

        KB_f_x = fns_KBvalue.fn_rms_kb(t, KBx_ifft*1000 , 0.125);
        KB_f_y = fns_KBvalue.fn_rms_kb(t, KBy_ifft*1000 , 0.125);
        KB_f_z = fns_KBvalue.fn_rms_kb(t, KBz_ifft*1000 , 0.125);

        max_KBx_f_mat(i_stn,:,i_flur)=max(KB_f_x);
        max_KBy_f_mat(i_stn,:,i_flur)=max(KB_f_y);
        max_KBz_f_mat(i_stn,:,i_flur)=max(KB_f_z);
        max_KB_f_mat(i_stn,:,i_flur) = max([max(KB_f_x); max(KB_f_y); max(KB_f_z)], [], 1);

        for ilb=1:num_lb

            % find the freq of max Vss
            Vss_Xvect=abs(Vss_Xmat(:,ilb));
            [~, index] = max(Vss_Xvect);
            f_x_max(i_stn,ilb,i_flur)=freq(index);
            Vss_Yvect=abs(Vss_Ymat(:,ilb));
            [~, index] = max(Vss_Yvect);
            f_y_max(i_stn,ilb,i_flur)=freq(index);
            Vss_Zvect=abs(Vss_Zmat(:,ilb));
            [~, index] = max(Vss_Zvect);
            f_z_max(i_stn,ilb,i_flur)=freq(index);


            %find the time of max Kb_f
            KB_f_x_vect=KB_f_x(:,ilb);
            [max_x, index] = max(KB_f_x_vect);
            t_x_max(i_stn,ilb,i_flur)=t(index);

            KB_f_y_vect=KB_f_y(:,ilb);
            [max_y, index] = max(KB_f_y_vect);
            t_y_max(i_stn,ilb,i_flur)=t(index);

            KB_f_z_vect=KB_f_z(:,ilb);
            [max_z, index] = max(KB_f_z_vect);
            t_z_max(i_stn,ilb,i_flur)=t(index);
            % Determine the maximum value among the three directions and select the corresponding time
            [~, max_idx] = max([max_x, max_y, max_z]);
            if max_idx == 1
                t_max(i_stn,ilb,i_flur) = t_x_max(i_stn,ilb,i_flur);
            elseif max_idx == 2
                t_max(i_stn,ilb,i_flur) = t_y_max(i_stn,ilb,i_flur);
            else
                t_max(i_stn,ilb,i_flur) = t_z_max(i_stn,ilb,i_flur);
            end

            % find the freq of max KBss
            KBss_Xvect=abs(KBss_Xmat(:,ilb));
            [max_x, ~] = max(KBss_Xvect);
            KBss_Yvect=abs(KBss_Ymat(:,ilb));
            [max_y, ~] = max(KBss_Yvect);
            KBss_Zvect=abs(KBss_Zmat(:,ilb));
            [max_z, ~] = max(KBss_Zvect);
            % Determine the maximum value among the three directions and
            % select the corresponding freq
            [~, max_idx] = max([max_x, max_y, max_z]);
            if max_idx == 1
                f_max(i_stn,ilb,i_flur) = f_x_max(i_stn,ilb,i_flur);
            elseif max_idx == 2
                f_max(i_stn,ilb,i_flur) = f_y_max(i_stn,ilb,i_flur);
            else
                f_max(i_stn,ilb,i_flur) = f_z_max(i_stn,ilb,i_flur);
            end
        end
    end
end
% [y_lim_x_D2,y_lim_y_D2,y_lim_z_D2,x_lim_D2,x_lim_x_D3,...
%     y_lim_x_D3,x_lim_y_D3,y_lim_y_D3,x_lim_z_D3,y_lim_z_D3]=...
%     fns_unitgeomdb.get_ylim(rf_fldr,name_evnt);
[y_lim_D2,x_lim_D2,x_lim_x_D3,...
    y_lim_x_D3,x_lim_y_D3,y_lim_y_D3,x_lim_z_D3,y_lim_z_D3]=...
    fns_unitgeomdb.get_ylim1(rf_fldr,name_evnt);
%%
ylbl_vect={'$\max[v_x(t)]$,~mm/s', '$\max[v_y(t)]$,~mm/s', '$\max[v_z(t)]$,~mm/s'};

ylbl=ylbl_vect{1};
cmp=cmpt{1};
fns_unitgeomdb.plot_DIN4150_3_XYZ1(v_ref,n_str+1,f_x_max,max_Vxmat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim_x_D3,y_lim_x_D3,rf_fldr)
ylbl=ylbl_vect{2};
cmp=cmpt{2};
fns_unitgeomdb.plot_DIN4150_3_XYZ1(v_ref,n_str+1,f_y_max,max_Vymat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim_y_D3,y_lim_y_D3,rf_fldr)
ylbl=ylbl_vect{3};
cmp=cmpt{3};
fns_unitgeomdb.plot_DIN4150_3_XYZ1(v_ref,n_str+1,f_z_max,max_Vzmat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim_z_D3,y_lim_z_D3,rf_fldr)

%%
bldtyp="ReinesWohngebiet";
ylbl_vect={'$\max[KB_{F}(x)]$', '$\max[KB_{F}(y)]$', '$\max[KB_{F}(z)]$'};
xlbl_vect={'$t_p$,~s','$f_d$,~Hz'};
% ylbl=ylbl_vect{1};
% xlbl=xlbl_vect{1};
% cmp=cmpt{1};
% fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_x_max,max_KBx_f_mat,V_s,n_stns,stn_vect,name_evnt,xlbl,ylbl,cmp,bldtyp,t_in,x_lim_D2,y_lim_x_D2,rf_fldr)
%
% ylbl=ylbl_vect{2};
% xlbl=xlbl_vect{1};
% cmp=cmpt{2};
% fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_y_max,max_KBy_f_mat,V_s,n_stns,stn_vect,name_evnt,xlbl,ylbl,cmp,bldtyp,t_in,x_lim_D2,y_lim_y_D2,rf_fldr)
%
% ylbl=ylbl_vect{3};
% xlbl=xlbl_vect{1};
% cmp=cmpt{3};
% fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_z_max,max_KBz_f_mat,V_s,n_stns,stn_vect,name_evnt,xlbl,ylbl,cmp,bldtyp,t_in,x_lim_D2,y_lim_z_D2,rf_fldr)
%%
% ylbl='$\max[KB_{F}]$';
% xlbl=xlbl_vect{2};
% cmp='max';
% 
% fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,f_max,max_KB_f_mat,V_s,n_stns,stn_vect,name_evnt,xlbl,ylbl,cmp,bldtyp,x_lim_D2,y_lim_D2,rf_fldr)
