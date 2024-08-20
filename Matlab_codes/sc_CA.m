%% Saves the results for the maximum values velocity and KB levels 
% at selected locations of the building models. The corresponding values
% dominant frequencies and time, respectively where the maximum v and KB
% are observed are also saved but they are not used later on. 
% (this remains to be modified)

%% Initialization
clear;clc;close all;

rf_fldr = fns_Inpt_BldPara.selectRfFldr();
disp(['selectedRfFldr: ', rf_fldr])
%----------------------------------------------------%
bld_soil_fndn = fns_Inpt_BldPara.select_bldsoilpara();
disp(['bld_soil_fndnPara: ', bld_soil_fndn])
%----------------------------------------------------%
data_set= fns_EvntData.select_event_stn();
disp(['selected_dataset: ', data_set])
% [stn_vect,date_evnt,time_evnt,r_vect]=fns_data_process.get_event_fordataprocess(data_set,evnt);
[l_vect,b_vect,h,wall_config,dampg_vect,bldcases]=fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);

[n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);

[event_vect,nzero,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=...
    fns_EvntData.getEventList(data_set);


%%
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};
%
% tic;
%% Start a parallel pool
% if isempty(gcp('nocreate'))
%     parpool;
% end

for ie = 1:length(event_vect)
    evnt=event_vect{ie};
    disp(['evnt: ', evnt])
    [stn_vect,date_evnt,time_evnt,r_vect]=fns_data_process.get_event_fordataprocess(data_set,evnt);
    %%
    n_stns=length(stn_vect);
    max_Vxmat = zeros(n_stns, bldcases, n_str + 1);
    max_Vymat = zeros(n_stns, bldcases, n_str + 1);
    max_Vzmat = zeros(n_stns, bldcases, n_str + 1);
    max_Vx_KB_f_mat = zeros(n_stns, bldcases, n_str + 1);
    max_Vy_KB_f_mat = zeros(n_stns, bldcases, n_str + 1);
    max_Vz_KB_f_mat = zeros(n_stns, bldcases, n_str + 1);
    f_x_max = zeros(n_stns, bldcases, n_str + 1);
    f_y_max = zeros(n_stns, bldcases, n_str + 1);
    f_z_max = zeros(n_stns, bldcases, n_str + 1);
    t_x_max = zeros(n_stns, bldcases, n_str + 1);
    t_y_max = zeros(n_stns, bldcases, n_str + 1);
    t_z_max = zeros(n_stns, bldcases, n_str + 1);
    for i_stn=1:n_stns
        stn=stn_vect{i_stn};
        disp(['stn: ', stn])
        [ff_fldr]=fns_EvntData.get_ff_fldr1(data_set,stn,date_evnt,time_evnt);
        %% Importing sensor data
        [f_inpt,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
            fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
            stn,date_evnt, time_evnt,n_snr,ff_fldr,cols);
        [f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
            fns_imprtdata.get_ff_inpt(bf_nm_v,s_dir,...
            stn,date_evnt, time_evnt,n_snr,ff_fldr,cols);
        % fns_plot.plt_ff(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
        %     stn,date_evnt, time_evnt,'Velocity~(m/s)','initial')

        [t_in,ff_Vt]=...
            fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
            stn,date_evnt, time_evnt,n_snr,ff_fldr,cols_t);
        t_in=t_in{1};

        n_c = length(cmpt);
        v_ref=5e-8;

        %% Importing Transfer function computed from ANSYS
        % Assuming 'n_str' is the number of stories
        Vss_zCell = cell(1, n_str + 1);
        Vss_xCell = cell(1, n_str + 1);
        Vss_yCell = cell(1, n_str + 1);
        V_KB_freq_zCell = cell(1, n_str + 1);
        V_KB_freq_xCell = cell(1, n_str + 1);
        V_KB_freq_yCell = cell(1, n_str + 1);

        for i_str = 0:n_str
            [f_vect, TF_amp_mat, TF_cpmlx_mat] = fns_Wall_and_DR.get_TF(rf_fldr,...
                n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f,...
                wall_config, dampg_vect, i_str, n_c);
            % Assuming 'n_c' is the number of components
            TFcpmlx_intrp = cell(1, n_c);
            Ucmplx_mat = cell(1, n_c);
            Vabs_mat = cell(1, n_c);
            Vcmplx_mat = cell(1, n_c);
            V_KB_freq = cell(1, n_c);

            for i_c = 1:n_c
                %% Calculating Velocity
                TFcpmlx_intrp{i_c}=interp1(f_vect,TF_cpmlx_mat{i_c},...
                    f_inpt{i_c},'linear','extrap');

                Ucmplx_mat{i_c}=ff_Ucmplx_mat{i_c}.*TFcpmlx_intrp{i_c};

                Vabs_mat{i_c}=abs(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
                Vcmplx_mat{i_c}=(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
                %% Apply filter transfer to KB(f)
                V_KB_freq{i_c} = fns_KBvalue.fn_kb_highpass(f_inpt{i_c},Vcmplx_mat{i_c},5.6);
            end
            % Original data
            Vss_zCell{i_str+1} = Vcmplx_mat{3};
            Vss_xCell{i_str+1} = Vcmplx_mat{1};
            Vss_yCell{i_str+1} = Vcmplx_mat{2};

            %Data with filter: KB(f)
            V_KB_freq_zCell{i_str+1} = V_KB_freq{3};
            V_KB_freq_xCell{i_str+1} = V_KB_freq{1};
            V_KB_freq_yCell{i_str+1} = V_KB_freq{2};
        end

        %% IFFT compute v(t)
        % Fs = 200;           %sampling rate
        Fs = 1./(t_in(3)-t_in(2));
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
        % Assuming 'n_str' for the number of stories
        Vz_ifft_cell = cell(1, n_str + 1);
        Vx_ifft_cell = cell(1, n_str + 1);
        Vy_ifft_cell = cell(1, n_str + 1);
        Vz_KB_ifft_cell = cell(1, n_str + 1);
        Vx_KB_ifft_cell = cell(1, n_str + 1);
        Vy_KB_ifft_cell = cell(1, n_str + 1);

        for i_flur = 1:n_str+1
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
            V_KB_Zmat = V_KB_freq_zCell{i_flur};
            V_KB_Xmat = V_KB_freq_xCell{i_flur};
            V_KB_Ymat = V_KB_freq_yCell{i_flur};
            Vz_KB_fft_pad = [V_KB_Zmat; conj(flipud(V_KB_Zmat(2:end-1,:)))];
            Vx_KB_fft_pad = [V_KB_Xmat; conj(flipud(V_KB_Xmat(2:end-1,:)))];
            Vy_KB_fft_pad = [V_KB_Ymat; conj(flipud(V_KB_Ymat(2:end-1,:)))];

            Vz_KB_ifft = ifft(Vz_KB_fft_pad*Fs, nfft, 1, 'symmetric');
            Vx_KB_ifft = ifft(Vx_KB_fft_pad*Fs, nfft, 1, 'symmetric');
            Vy_KB_ifft = ifft(Vy_KB_fft_pad*Fs, nfft, 1, 'symmetric');
            Vz_KB_ifft = Vz_KB_ifft(1:length(t_in),:);
            Vx_KB_ifft = Vx_KB_ifft(1:length(t_in),:);
            Vy_KB_ifft = Vy_KB_ifft(1:length(t_in),:);

            % Store the result in Cell : {floor num}
            Vz_KB_ifft_cell{i_flur} = Vz_KB_ifft;
            Vx_KB_ifft_cell{i_flur} = Vx_KB_ifft;
            Vy_KB_ifft_cell{i_flur} = Vy_KB_ifft;

            t = (0:length(Vz_ifft)-1) / Fs;

            max_Vxmat(i_stn,:,i_flur)=max(Vx_ifft);
            max_Vymat(i_stn,:,i_flur)=max(Vy_ifft);
            max_Vzmat(i_stn,:,i_flur)=max(Vz_ifft);

            KB_f_x = fns_KBvalue.fn_rms_kb(t, Vx_KB_ifft*1000 , 0.125);
            KB_f_y = fns_KBvalue.fn_rms_kb(t, Vy_KB_ifft*1000 , 0.125);
            KB_f_z = fns_KBvalue.fn_rms_kb(t, Vz_KB_ifft*1000 , 0.125);

            max_Vx_KB_f_mat(i_stn,:,i_flur)=max(KB_f_x);
            max_Vy_KB_f_mat(i_stn,:,i_flur)=max(KB_f_y);
            max_Vz_KB_f_mat(i_stn,:,i_flur)=max(KB_f_z);

            for ilb=1:bldcases

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
                [~, index] = max(KB_f_x_vect);
                t_x_max(i_stn,ilb,i_flur)=t(index);

                KB_f_y_vect=KB_f_y(:,ilb);
                [~, index] = max(KB_f_y_vect);
                t_y_max(i_stn,ilb,i_flur)=t(index);

                KB_f_z_vect=KB_f_z(:,ilb);
                [~, index] = max(KB_f_z_vect);
                t_z_max(i_stn,ilb,i_flur)=t(index);
            end
        end
    end
    fns_CloudAnalysis.save_DINvals(data_set,date_evnt,time_evnt,...
        rf_fldr,bld_soil_fndn,f_x_max,f_y_max,f_z_max,...
        max_Vxmat,max_Vymat,max_Vzmat,t_x_max,t_y_max,t_z_max,...
        max_Vx_KB_f_mat,max_Vy_KB_f_mat,max_Vz_KB_f_mat)


end
% t_ex=toc;
% disp(['Cloud analysis: evnt processing time=',t_ex])