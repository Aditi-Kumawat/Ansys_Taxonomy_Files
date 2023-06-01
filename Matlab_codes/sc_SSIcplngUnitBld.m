%%
clear;clc;
close all;
% Define the number of storeys, rooms in x-y-direction
n_str = 1;
n_rx = 1;
n_ry = 1;
% Define the length, width, and height of the building
l = 4;
b = 4;
h = 3;
% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'FOOTING';
V_s =40;
rho_s=1260;
% Define the size of the elements
n_esize = 0.25;
if strcmp(ftyp,'PLATE')
    B_f = n_esize/2;
    L_f = n_esize/2;
else
    B_f = 0.75;
    L_f = 0.75;
end
t_f=0.2;
n_f=4;
vol_fndn=4*L_f*B_f*t_f*n_f;
rho_bld=2500;

vS_vect =[40 155 270 385 500];
rhos_vect=[1260,1570,1880,2190,2500];
%% Importing Transfer Function
rf_fldr = 'Results_SSIcplng';
bf_nm = 'Disp_Center_%s_%d_l%d_b%d';
cols1 = {'Freq', 'AMPL','PHASE','REAL','IMAG'};
cmpt = {'X', 'Y', 'Z'};
n_c = length(cmpt);

%%
%floor 0 values are all 1 as the load is applied at the footing
i_str = 1;
[f_out,Uamp_mat,Ucpmlx_mat]=fns_imprtdata.get_TF(...
    n_str, n_rx, n_ry,...
    l, b, ftyp, V_s, L_f, B_f,...
    bf_nm,i_str,cmpt,n_c,rf_fldr,cols1);

%% Recheck the parameters for K;C,and M
k_ratio=1; %Ratio of B/C
BLratio_vect=[1 2 4];
BLratio=BLratio_vect(k_ratio);
mfXY=0.5*vol_fndn*rho_bld;
nu=1/3;

G_s=rho_s.*V_s.^2;
KY=((2*G_s*L_f)/(2-nu))*(2+2.5*(BLratio^0.85));
cY=((L_f./V_s).*0.58.*KY);
KZ=((2*G_s*L_f)/(1-nu))*(0.73+1.54*(BLratio^0.75));
cZ=((L_f./V_s).*0.85.*KZ);

G_s_vect=rhos_vect.*vS_vect.^2;
KX_vect=((2*G_s_vect*L_f)/(2-nu))*(2+2.5*(BLratio^0.85));
cX_vect=((L_f./vS_vect).*0.58.*KX_vect);
KY_vect=((2*G_s_vect*L_f)/(2-nu))*(2+2.5*(BLratio^0.85));
cY_vect=((L_f./vS_vect).*0.58.*KY_vect);
KZ_vect=((2*G_s_vect*L_f)/(1-nu))*(0.73+1.54*(BLratio^0.75));
cZ_vect=((L_f./vS_vect).*0.85.*KZ_vect);
KF_mat=[KX_vect;KY_vect;KZ_vect];
cF_mat=[cX_vect;cY_vect;cZ_vect];
mf=mfXY;
%%
f_vect1=f_out;
u_Rmat=zeros(length(vS_vect),length(f_vect1));
for i_c=1:3
    U_bld=Ucpmlx_mat{i_c};
    for i_v=1:length(vS_vect)
        K1=KF_mat(i_c,i_v);
        c1=cF_mat(i_c,i_v);
        omg=2*pi*f_vect1;
        u1_vect=1./(-mf*omg.^2+1i*omg.*c1+K1);
        u_R=(u1_vect)./(U_bld);
        u_Rmat(i_v,:)=abs(u_R);
    end

    figure
    hold on
    leg_vect = cell(1, length(vS_vect)-1);
    for i_v = 1:length(vS_vect)
        plot(f_vect1, (u_Rmat(i_v,:)), 'LineWidth', 1.5)
        leg_vect{i_v} = sprintf('$V_s =$ %d m/s', vS_vect(i_v));
    end

    legend(leg_vect, 'Box', 'off', 'Interpreter', 'latex',...
        'FontSize', 12)
    xlabel({'Frequency (Hz)'}, 'FontSize', 12,...
        'Interpreter', 'latex')
    ylabel('$K_{Bld}/K_{Foundation}$', 'FontSize', 12,...
        'Interpreter', 'latex')
    txt = ['Component:~',num2str(i_c)];
    title(txt,'Interpreter','latex','FontSize',11);
    set(gca, 'XTickLabelMode', 'auto');
    set(gca, 'YTickLabelMode', 'auto');
    box on
    set(gcf, 'Units', 'inches', 'Position', [30 4 6 3],...
        'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125]);
    filename = ['URatio_UnitBld_abs_', num2str(i_c),...
        '_l', num2str(l), '_by_b', num2str(b),...
        '_ftyp_', ftyp, '_Vsmin_', num2str(V_s),...
        '_Lf_', num2str(L_f), '_Bf_', num2str(B_f), '.png'];

    cd SAVE_FIGS
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end
    saveas(gcf, fullfile(rf_fldr, filename));
    cd ..
    cd ..
    cd Matlab_codes
end



