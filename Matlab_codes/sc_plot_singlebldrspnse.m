%%
clear;clc;
close all;
% Define the number of storeys, rooms in x-y-direction
n_str = 3;
n_rx = 2;
n_ry = 3;
% Define the length, width, and height of the building
l = 5;
b = 5;
h = 3;
% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'PLATE';
% Define the velocity of the excitation
V_s =450;
% Define the size of the elements
n_esize = 0.5;
if strcmp(ftyp,'PLATE')
    B_f = n_esize/2;
    L_f = n_esize/2;
else
    B_f = 0.75;
    L_f = 0.75;
end
%%
%%
bf_nm_vt = 'v_%d_%s_%s_%s';
cols_t = {'tim', 'val'};
%----------------------------------------------------%
name_evnt= fns_EvntData.select_event_stn();
disp(['selected_event: ', name_evnt])
[evnt,stn,~,date,time,nzero,ff_fldr,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=...
    fns_EvntData.get_event_stn(name_evnt);
[f_inpt,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
    fns_imprtdata.get_ff_inpt(bf_nm_u,s_dir,...
    stn,date, time,n_snr,ff_fldr,cols);
[t_in,ff_Vt]=...
    fns_imprtdata.get_ff_tim(bf_nm_vt,s_dir,...
    stn,date, time,n_snr,ff_fldr,cols_t);
t_in=t_in{1};
%%
rf_fldr = 'MultiUnitBld_GeomVary_3lby2';
% rf_fldr = 'UnitBld_GeomVary';
% rf_fldr = 'TF_F_Indpn';

bf_nm = 'Disp_Center_%s_%d_l%d_b%d';

cols = {'Freq', 'AMPL','PHASE','REAL','IMAG'};

cmpt = {'X', 'Y', 'Z'};

f_vect = [];

folder = fns_plot.get_fldrnm(n_str,n_rx,n_ry,l,b,ftyp,V_s,L_f,B_f);


%% IFFT compute v(t)
Fs = 200;           %sampling rate
nfft = length(t_in);
%frequnecy array
freq = Fs / 2 * linspace(0, 1, nfft/2+1); %length = nfft/2+1
dfz=freq(3)-freq(2);
%%
u_ref=1;
fct=1e3;
for i_str = 0:n_str
    fil_nm = arrayfun(@(x) sprintf(bf_nm, x{1}, i_str, l, b),...
        cmpt, 'UniformOutput', false);

    cd ..
    cd APDL_codes
    cd Results_Ansys
    fil_pth = fullfile(rf_fldr, folder, fil_nm);
    U_all = cellfun(@(x) readtable(x),fil_pth,'UniformOutput',false);

    for i_component = 1:3
        U = U_all{i_component};
        U.Properties.VariableNames = cols;
        f_vect = U.Freq;
        TFamp_mat{i_component} = U.AMPL;
        TFr_mat{i_component} = U.REAL;
        TFim_mat{i_component} = U.IMAG;
        TFcmplx_mat{i_component} = U.REAL+1i.*U.IMAG;
    end
    cd ..
    cd ..
    cd Matlab_codes
    for i_c = 1:3
        %% Calculating Velocity
        TFcpmlx_intrp{i_c}=interp1(f_vect,TFcmplx_mat{i_c},...
            f_inpt{i_c},'linear','extrap');
        Ucmplx_mat{i_c}=ff_Ucmplx_mat{i_c}.*TFcpmlx_intrp{i_c};

        Vabs_mat{i_c}=abs(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
        Vcmplx_mat{i_c}=(Ucmplx_mat{i_c}.*1i.*f_inpt{i_c}*2*pi);
    end
    Vamp_X = Vabs_mat{1};
    Vamp_Y = Vabs_mat{2};
    Vamp_Z = Vabs_mat{3};
    % Original data
    Vss_zCell{i_str+1} = Vcmplx_mat{3};
    Vss_xCell{i_str+1} = Vcmplx_mat{1};
    Vss_yCell{i_str+1} = Vcmplx_mat{2};
    %% Plotting
    ha_cl = @colors;
    lStyl = {'-', '--', ':', '-.'};
    lcol = {ha_cl('boston university red'),ha_cl('black'),...
        ha_cl('denim')};

    figure
    plot(f_inpt{1}, Vamp_X*fct, 'LineStyle', lStyl{1}, 'Color', lcol{1},...
        'DisplayName', 'X-dir', 'LineWidth', 1)
    hold on
    plot(f_inpt{1}, Vamp_Y*fct, 'LineStyle', lStyl{1}, 'Color', lcol{2},...
        'DisplayName', 'Y-dir', 'LineWidth', 1)
    hold on
    plot(f_inpt{1}, Vamp_Z*fct, 'LineStyle', lStyl{4}, 'Color', lcol{3},...
        'DisplayName', 'Z-dir', 'LineWidth', 0.5)

    legend('show', 'Box', 'off', 'Interpreter', 'latex',...
        'FontSize', 8)
    xlabel({'f, Hz'}, 'FontSize', 10,...
        'Interpreter', 'latex')
    ylabel('$v$, mm/s', 'FontSize', 10,...
        'Interpreter', 'latex')

    set(gca, 'XTickLabelMode', 'auto');
    set(gca, 'YTickLabelMode', 'auto');
    set(gca,'FontSize',8, 'Box', 'on','LineWidth',0.2,...
        'TickLabelInterpreter','latex',...
        'TickLength',[0.01,0.01]);
    set(gcf, 'Units', 'inches', 'Position',...
        [18 3 3.0 4/3], 'PaperUnits', 'Inches',...
        'PaperSize', [3.0 4/3]);
    % if i_str==0
    %     ylim([0.5,1.6])
    %     % else
    %     %     ylim([0,15.3])
    % end
    xlim([0,80])
    filename = ['VcentXYZ_', num2str(i_str),...
        '_n_rooms_X_', num2str(n_rx),...
        '_n_rooms_Y_', num2str(n_ry),...
        '_l', num2str(l), '_by_b', num2str(b),...
        '_ftyp_', ftyp, '_Vs_', num2str(V_s),...
        '_Lf_', num2str(L_f), '_Bf_', num2str(B_f),stn, '.pdf'];

    cd SAVE_FIGS
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end
    saveas(gcf, fullfile(rf_fldr, filename));
    cd ..
    cd ..
    cd Matlab_codes
end

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
    t = (0:length(Vz_ifft)-1) / Fs;

    %% Plotting
    ha_cl = @colors;
    lStyl = {'-', '--', ':', '-.'};
    lcol = {ha_cl('boston university red'),ha_cl('black'),...
        ha_cl('denim')};

    figure
    plot(t, Vx_ifft*fct, 'LineStyle', lStyl{1}, 'Color', lcol{1},...
        'DisplayName', 'X-dir', 'LineWidth', 1)
    hold on
    plot(t, Vy_ifft*fct, 'LineStyle', lStyl{1}, 'Color', lcol{2},...
        'DisplayName', 'Y-dir', 'LineWidth', 1)
    hold on
    plot(t, Vz_ifft*fct, 'LineStyle', lStyl{4}, 'Color', lcol{3},...
        'DisplayName', 'Z-dir', 'LineWidth', 0.5)

    legend('show', 'Box', 'off', 'Interpreter', 'latex',...
        'FontSize', 8)
    xlabel({'$t$, s'}, 'FontSize', 10,...
        'Interpreter', 'latex')
    ylabel('$v$, mm/s', 'FontSize', 10,...
        'Interpreter', 'latex')

    set(gca, 'XTickLabelMode', 'auto');
    set(gca, 'YTickLabelMode', 'auto');
    set(gca,'FontSize',8, 'Box', 'on','LineWidth',0.2,...
        'TickLabelInterpreter','latex',...
        'TickLength',[0.01,0.01]);
    set(gcf, 'Units', 'inches', 'Position',...
        [18 3 3.0 4/3], 'PaperUnits', 'Inches',...
        'PaperSize', [3.0 4/3]);
    % if i_str==0
    %     ylim([0.5,1.6])
    %     % else
    %     %     ylim([0,15.3])
    % end
    xlim([0,40])
    filename = ['Vt_', num2str(i_flur),...
        '_n_rooms_X_', num2str(n_rx),...
        '_n_rooms_Y_', num2str(n_ry),...
        '_l', num2str(l), '_by_b', num2str(b),...
        '_ftyp_', ftyp, '_Vs_', num2str(V_s),...
        '_Lf_', num2str(L_f), '_Bf_', num2str(B_f),stn, '.pdf'];

    cd SAVE_FIGS
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end
    saveas(gcf, fullfile(rf_fldr, filename));
    cd ..
    cd ..
    cd Matlab_codes

end