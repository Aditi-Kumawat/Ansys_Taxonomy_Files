
%% Initialization
clear;clc;close all;

% Define the number of storeys, rooms in x-,y-direction
n_str = 2;
n_rx = 2;
n_ry = 3;

h=3;
l=5;
b=3;
% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'PLATE';
% Define the velocity of the excitation
V_s = 450;
n_esize_vect =0.5;

% Define the size of the elements
%% Importing Transfer Function
% Define the name of the folder where the results are stored
rf_fldr = 'nele_Vary';
bf_nm = 'Disp_Center_%s_%d_nel%d';

% Define the column names in the results table
cols1 = {'Freq', 'AMPL','PHASE','REAL','IMAG'};

% Define the components to extract from the results
cmpt = {'X', 'Y', 'Z'};

%% Importing Data
rec=1324;

ff_fldr = 'GM/SeisSol/data_noVs30_Bld_new';
% Define the base file name for the U center results in ANSYS
bf_nm_u = 'fftd_%d_%s_%d';
bf_nm_v = 'fftv_%d_%s_%d';
stn='rec';
% Define the column names in the results table
cols = {'Freq', 'Re','Im','Amp'};

% Define the components to extract from the results
s_dir = [1 2 3];
% length of free field data for each sensor is different
n_snr = numel(s_dir);

[f_inpt_U,ff_Uamp_mat,ff_Ur_mat,ff_UIm_mat,ff_Ucmplx_mat]=...
    fns_imprtdata.get_inpt_seisl(bf_nm_u,s_dir,...
    stn,rec,n_snr,ff_fldr,cols);
[f_inpt_V,ff_Vamp_mat,ff_Vr_mat,ff_VIm_mat,ff_Vcmplx_mat]=...
    fns_imprtdata.get_inpt_seisl(bf_nm_v,s_dir,...
    stn,rec,n_snr,ff_fldr,cols);
%%
ss_fldr = 'GM/SeisSol/data_noVs30_Bld_new';
rec_b_vect=[1324 1325 1326];
[bldUampmat,bldUcmplx_mat]=fns_imprtdata.get_bldata_seisl(bf_nm_u,...
    stn,rec_b_vect,s_dir,ss_fldr,cols);
[bldVampmat,bldVcmplx_mat]=fns_imprtdata.get_bldata_seisl(bf_nm_v,...
    stn,rec_b_vect,s_dir,ss_fldr,cols);

% Plotting Data
date='1324';
time='1324';
fns_plot.plt_ff(f_inpt_V, ff_Vamp_mat,bf_nm_v,123,...
    stn,date, time,'Velocity~(m/s)','initial')
% fns_plot.plt_ff(f_inpt_U, ff_Uamp_mat,bf_nm_u,123,...
%     stn,date, time,'Displacement~(m)','initial')

%%
n_c = length(cmpt);
vel_abs_mat_Z_cell = cell(1, n_str+1);
fref_cmplx_mat_1=cell(1, n_c);
DISP_cmplx_mat=cell(1, n_c);
VEL_abs_mat=cell(1, n_c);
fref_vel_amp_mat_1=cell(1, n_c);
Ucmplx_mat=cell(1, n_c);
Vabs_mat=cell(1, n_c);
Vcmplx_mat=cell(1, n_c);
for i_str = 0:n_str
    bld=bldVampmat{i_str+1};
    for i_n=1:length(n_esize_vect)
        n_esize= n_esize_vect(i_n);
        n_esize_name=n_esize*1000;
        if strcmp(ftyp,'PLATE')
            B_f = n_esize/2;
            L_f = n_esize/2;
        else
            B_f = 0.75;
            L_f = 0.75;
        end
        fldr = fns_plot.get_fldrnm(n_str, n_rx, n_ry,...
            l, b,ftyp, V_s, L_f, B_f);

        fil_nm = arrayfun(@(x) sprintf(bf_nm, x{1}, i_str, n_esize_name),...
            cmpt, 'UniformOutput', false);

        cd ..
        cd Results_Ansys
        fil_pth = fullfile(rf_fldr, fldr, fil_nm);
        U_all = cellfun(@(x) readtable(x), fil_pth,...
            'UniformOutput', false);
        cd ..
        cd Matlab_codes

        for i_c = 1:n_c
            Uc = U_all{i_c};
            Uc.Properties.VariableNames = cols1;
            if i_n == 1
                f_vect = Uc.Freq;
                TFamp{i_c} = Uc.AMPL;
                UR_mat{i_c} = Uc.REAL;
                UIm_mat{i_c} = Uc.IMAG;
                Ucpmlx_mat{i_c} = Uc.REAL+1i.*Uc.IMAG;
            else
                TFamp{i_c}(:, i_n)= Uc.AMPL;
                UR_mat{i_c}(:, i_n)= Uc.REAL;
                UIm_mat{i_c}(:, i_n)= Uc.IMAG;
                Ucpmlx_mat{i_c}(:, i_n)=Uc.REAL+1i.*Uc.IMAG;
            end
            %% Calculating Velocity
            Ucpmlx_intrp{i_c}=interp1(f_vect,Ucpmlx_mat{i_c},...
                f_inpt_U{i_c},'linear','extrap');
            Ucmplx_mat{i_c}=ff_Ucmplx_mat{i_c}.*Ucpmlx_intrp{i_c};

            Vabs_mat{i_c}=abs(Ucmplx_mat{i_c}.*1i.*f_inpt_U{i_c}*2*pi);
            Vcmplx_mat{i_c}=(Ucmplx_mat{i_c}.*1i.*f_inpt_U{i_c}*2*pi);

        end
    end
    for i_c = 1:n_c
        %% Plotting Transfer Function
        ha_col = @colors;
        lStyl = {'-', ':', '--', '-.'};
        lcol = {ha_col('boston university red'), ha_col('cadmium green'),...
            ha_col('black'),ha_col('denim'),...
            ha_col('dark goldenrod')};
        ustr_vect={'$u_x$~(m)','$u_y$~(m)','$u_z$~(m)'};
        figure
        for i_n1=1:length(n_esize_vect)
            n_esize= n_esize_vect(i_n1);
            txt_1 = ['element~size:~',num2str(n_esize),'m'];
            if i_c==1
                i_col=1;
            elseif i_c==2
                i_col=2;
            else
                i_col=4;
            end
            hold on
            plot(f_inpt_U{i_c},Vabs_mat{i_c}(:,i_n1),...
                'linestyle',lStyl{mod(i_n1-1,numel(lStyl))+1},...
                'DisplayName',txt_1,'LineWidth',1.5,...
                'Color',lcol{mod(i_n1-1,numel(lcol))+1})

        end
        hold on
        plot(f_inpt_U{i_c}, bld(:,i_c),'DisplayName','SeisSol',...
            'Color','k', 'LineWidth', 1.2);


        fns_plot.setPltProps(ustr_vect,i_c);

        filename = ['TFcmpr_walls',cmpt{i_c},...
            num2str(i_str),'_n_rooms_X_', num2str(n_rx),...
            '_n_rooms_Y_', num2str(n_ry),...
            '_ftyp_', ftyp,'_Vs_', num2str(V_s),...
            '_Lf_', num2str(L_f), '_Bf_', num2str(B_f),'.pdf'];

        %         cd SAVE_FIGS
        %         if ~exist(rf_fldr_wall, 'dir')
        %             mkdir(rf_fldr_wall);
        %         end
        % %         saveas(gcf, fullfile(rf_fldr_wall, filename));
        %         cd ..
        %         cd ..
        %         cd Matlab_codes
    end
end
