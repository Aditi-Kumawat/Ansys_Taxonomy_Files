%% Initialization
clear;clc;close all;

% Define the number of storeys, rooms in x- y-direction
n_str = 3;
n_rx = 2;
n_ry = 3;

% Define the length, width, and height of the building
l_vect=5;
b_vect=5;
h = 3;
%
% wall_config = [1,2,3,4,5,6,7,8,9,10];
dampg_vect = 0:0.005:0.1;
num_dampg=length(dampg_vect);

%
% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'PLATE';

% Define the velocity of the excitation
V_s = 450;

% Define the size of the elements
n_esize = 0.5;

% Calculate the length and width of the footing based on the
% foundation type
if strcmp(ftyp,'PLATE')
    B_f = n_esize/2;
    L_f = n_esize/2;
else
    B_f = 0.75;
    L_f = 0.75;
end
%%
name_evnt='Poing';
%%
if strcmp(name_evnt, 'Poing')
    evnt='Po2016';
    stn_vect={'POI01', 'POI02', 'POI03'};
    %stn_vect={'POI01'};
    date='2016_12_20';
    time='03_30_51';
elseif strcmp(name_evnt, 'Unterhaching')
    evnt='Part1';
    stn_vect={'UH1', 'UH2', 'UH3'};
    date='2013_04_16';
    time='21_51_42';
end

n_stns=length(stn_vect);

%% Importing Data
bf_nm_u = 'fftd_%d_%s_%s_%s';
bf_nm_v = 'fftv_%d_%s_%s_%s'
nzero=1;   % what is nzero for? nzero = 6
cols = {'Freq', 'Re','Im','Amp'};
s_dir = [1 2 3];
n_snr = numel(s_dir);
%% Importing Transfer Function
rf_fldr = 'Vary_DampRatio';
bf_nm = 'Disp_Center_%s_%d_l%d_b%d';
cols1 = {'Freq', 'AMPL','PHASE','REAL','IMAG'};
cmpt = {'X', 'Y', 'Z'};

%%
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};

for i_stn=1:n_stns
    stn=stn_vect{i_stn};
    display(stn)
    if strcmp(name_evnt, 'Poing')
        ff_fldr = fullfile('GM','GM_POI2016',stn);
    elseif strcmp(name_evnt, 'Unterhaching')
        fldr_nm = [stn, '_', evnt];
        ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
    end
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
    Vabs_zCell = cell(1, n_str+1);
    fref_cmplx_mat_1=cell(1, n_c);
    DISP_cmplx_mat=cell(1, n_c);
    VEL_abs_mat=cell(1, n_c);
    fref_vel_amp_mat_1=cell(1, n_c);
    v_ref=5e-8;

    %% Importing Transfer function computed from ANSYS
    for i_str = 0:n_str
        [f_vect,TF_amp_mat,TF_cpmlx_mat]=...
            fns_Wall_and_DR.get_TF_DR(n_str, n_rx, n_ry,...
            l_vect, b_vect, ftyp, V_s, L_f, B_f,dampg_vect,...
            bf_nm,i_str,cmpt,n_c,rf_fldr,cols1);
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

        %Corresponding Frequency
        freq_=f_inpt_V{3}(nzero:end);
    end

    %% IFFT compute v(t)
    Fs = 200;           %sampling rate
    nfft = length(t_in);
    %frequnecy array
    freq = Fs / 2 * linspace(0, 1, nfft/2+1); %length = nfft/2+1
    dfz=freq(3)-freq(2);
    %% ----------------------------------------------------------------
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

        for ilb=1:num_dampg

            % find the freq of max Vss
            Vss_Xvect=abs(Vss_Xmat(:,ilb));
            [value, index] = max(Vss_Xvect);
            f_x_max(i_stn,ilb,i_flur)=freq(index);
            Vss_Yvect=abs(Vss_Ymat(:,ilb));
            [value, index] = max(Vss_Yvect);
            f_y_max(i_stn,ilb,i_flur)=freq(index);
            Vss_Zvect=abs(Vss_Zmat(:,ilb));
            [value, index] = max(Vss_Zvect);
            f_z_max(i_stn,ilb,i_flur)=freq(index);


            %find the time of max Kb_f
            KB_f_x_vect=KB_f_x(:,ilb);
            [value, index] = max(KB_f_x_vect);
            t_x_max(i_stn,ilb,i_flur)=t(index);

            KB_f_y_vect=KB_f_y(:,ilb);
            [value, index] = max(KB_f_y_vect);
            t_y_max(i_stn,ilb,i_flur)=t(index);

            KB_f_z_vect=KB_f_z(:,ilb);
            [value, index] = max(KB_f_z_vect);
            t_z_max(i_stn,ilb,i_flur)=t(index);
        end
    end
end

strtyp='DRVary';
ylbl_vect={'$v_{x,max}$,~m/s', '$v_{y,max}$,~m/s', '$v_{z,max}$,~m/s'};
x_lim=[0 20];
ylbl=ylbl_vect{1}
cmp=cmpt{1};
fns_unitgeomdb.plot_DIN4150_3_XYZ(v_ref,n_str+1,f_x_max,max_Vxmat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim,strtyp,rf_fldr)
x_lim=[0 20];
ylbl=ylbl_vect{2}
cmp=cmpt{2};
fns_unitgeomdb.plot_DIN4150_3_XYZ(v_ref,n_str+1,f_y_max,max_Vymat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim,strtyp,rf_fldr)
x_lim=[0 60];
ylbl=ylbl_vect{3}
cmp=cmpt{3};
fns_unitgeomdb.plot_DIN4150_3_XYZ(v_ref,n_str+1,f_z_max,max_Vzmat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,x_lim,strtyp,rf_fldr)

%%
bldtyp="ReinesWohngebiet";

ylbl_vect={'$KB_{Fmax}(x)$', '$KB_{Fmax}(y)$', '$KB_{Fmax}(z)$'};
ylbl=ylbl_vect{1}
cmp=cmpt{1};
y_lim=[0 1];
fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_x_max,max_Vx_KB_f_mat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,bldtyp,t_in,y_lim,strtyp,rf_fldr)

ylbl=ylbl_vect{2}
cmp=cmpt{2};
y_lim=[0 1];
fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_y_max,max_Vy_KB_f_mat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,bldtyp,t_in,y_lim,strtyp,rf_fldr)

ylbl=ylbl_vect{3}
cmp=cmpt{3};
y_lim=[0 1];
fns_unitgeomdb.plot_DIN4150_2_XYZ(n_str+1,t_z_max,max_Vz_KB_f_mat,V_s,n_stns,stn_vect,name_evnt,ylbl,cmp,bldtyp,t_in,y_lim,strtyp,rf_fldr)
