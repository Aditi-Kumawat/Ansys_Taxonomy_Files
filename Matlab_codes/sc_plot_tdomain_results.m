%% Initialization

clear;clc;
close all
% Define the number of storeys, rooms in x-y-direction

n_str = 2;
n_rx = 2;
n_ry = 3;

% Define the length, width, and height of the building
l_vect=5;
b_vect=3;
h = 3;

% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'PLATE';
% Define the foundation behaviour and analysis type
st_dtls='SSnoVs30Bld_DR_0pt05';
% st_dtls='SSnoVs30Bld_DR1002DR2004';
st_dtls='SS_HomogeneousBldDR1002DR2004_at';

% Define the velocity of the excitation
V_s = 1500;

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
% Define the name of the folder where the results are stored
rf_fldr = 'tDmn_input';

% Define the base file name for the U center results in ANSYS
bf_nm = 'Disp_Center_%s_%d_l%d_b%d';

% Define the column names in the results table
cols = {'time', 'val'};

% Define the components to extract from the results
cmpt = {'X', 'Y', 'Z'};
n_c = length(cmpt);

%% Importing and plotting Data for free-field vibrations in seissol (used as input for Ansys)
ff_fldr = 'GM/SeisSol/homo_soil/attn';
% Define the file name and path
fil_nm = 'vt_rec_17.txt';
% funs_plot_properties.import_tim_data(fref_data_folder,filename)
cd ..
% Combine the directory path and file name using the file separator
fil_pth = [ff_fldr filesep fil_nm];

% Import the data
data_seissol_freefield = readtable(fil_pth);
cd Matlab_codes

% Assign column names
data_seissol_freefield.Properties.VariableNames={'t_vect', 'data_x', 'data_y', 'data_z'};
%%

figure
for i_c=1:n_c
    subplot(n_c,1,i_c)
    plot(data_seissol_freefield.t_vect, data_seissol_freefield{:,i_c+1},'DisplayName',cmpt{i_c},...
        'LineWidth', 1.2)
    xlabel('Time (s)','Interpreter','latex','FontSize',12);
    ylabel('Velocity (m/s)','Interpreter','latex','FontSize',12);
    legend show
    legend('Box','off','Interpreter','latex','FontSize',11)
    %     title(sprintf('Component %d', i_c),'Interpreter','latex',...
    %         'FontSize',11);
    set(gcf,'Units','inches', 'Position', [18 3 4.5 6],...
        'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125]);
    fil_nm = ['FF_vel_t', '.png'];
    cd SAVE_FIGS
    % Create the directory (if it doesn't already exist)
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end

    % Save the figure in the results folder
    saveas(gcf, fullfile(rf_fldr, fil_nm));
    cd ..
    cd ..
    cd Matlab_codes
end
%% Importing data for building vibrations generated using Seisol
rec_vect=[18 19 20];



%% Importing the time-domain results generated in Ansys for disp and convert those to velocity
ut_mat=cell(1, n_c);
u_fft=cell(1, n_c);
vt_mat=cell(1, n_c);

for i_str = 0:n_str
    rec_f=rec_vect(i_str+1);
    fil_nm1 = strcat(sprintf('vt_rec_%d',rec_f), '.txt');
    % funs_plot_properties.import_tim_data(fref_data_folder,filename)
    cd ..
    % Combine the directory path and file name using the file separator
    fil_pth = [ff_fldr filesep fil_nm1];

    % Import the data
    data_seissol = readtable(fil_pth);
    cd Matlab_codes

    % Assign column names
    data_seissol.Properties.VariableNames = {'t', 'data_x', 'data_y', 'data_z'};
    for i_l=1:length(l_vect)
        l=l_vect(i_l);
        for i_b=1:length(b_vect)
            b=b_vect(i_b);
            if b>l
                break
            end
            fldr = fns_plot.get_fldrnm_rec(n_str, n_rx, n_ry,...
                l, b,ftyp, V_s, L_f, B_f,st_dtls);
            filNms = arrayfun(@(x) sprintf(bf_nm, x{1}, i_str, l, b),...
                cmpt, 'UniformOutput', false);

            cd ..
            cd APDL_codes
            cd Results_Ansys
            fil_pths = fullfile(rf_fldr, fldr, filNms);
            U_all = cellfun(@(x) readtable(x), fil_pths,...
                'UniformOutput', false);
            cd ..
            cd ..
            cd Matlab_codes

            for i_c = 1:n_c
                U = U_all{i_c};
                U.Properties.VariableNames = cols;
                %                 if i_l == 1
                t_vect = U.time;
                dt=t_vect(3)-t_vect(2);
                ut_mat{i_c} = U.val;
                %%
                vt_mat{i_c}=[0;diff(ut_mat{i_c})]./dt;

            end
        end
    end

    %% Plotting the comparison in time-domain
    sig_mat=vt_mat;
    figure
    for i_s = 1:n_c

        for i_col = 1:size(sig_mat{i_s}, 2)
            sig_vect = sig_mat{i_s}(:, i_col);
            % Plot the time-domain signal
            %         figure
            l=l_vect(i_col);
            b=b_vect(i_col);
            txt = ['SeisSol';'LPM-FEM'];
            %             txt = 'LPM-FEM';
            txt_2 = ['Floor:',num2str(i_str),...
                ',~Component:',num2str(i_s)];
            subplot(n_c, size(sig_mat{i_s}, 2),...
                (i_s-1)*size(sig_mat{i_s}, 2) + i_col);
            plot(data_seissol.t, data_seissol{:,i_s+1}, 'LineWidth', 1.5)
            hold on
            plot(t_vect-0.005, real(sig_vect),'LineStyle',':','LineWidth', 1.5);
            xlabel('Time (s)','Interpreter','latex','FontSize',12);
            ylabel('Velocity (m/s)','Interpreter','latex','FontSize',12);
            if i_s==1 && i_str==2
                ylim([-2.5e-4, 2.5e-4])
            elseif i_s==2 && i_str==2
                ylim([-2e-3, 2e-3])
            elseif i_s==3 && i_str==2
                ylim([-1e-3, 1e-3])
            end
            title(txt_2,'Interpreter','latex','FontSize',11);
            legend(txt);
            legend('Box','off','Interpreter','latex','FontSize',11)
            set(gcf,'Units','inches', 'Position', [18 3 4.5 6],...
                'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125]);
        end
    end
    fil_nm = ['Tdomain_vel_t_','_Floor_', num2str(i_str),...
        '_ftyp_', ftyp,'_Vs_', num2str(V_s),...
        '_Lf_', num2str(L_f), '_Bf_', num2str(B_f), '.emf'];
    cd SAVE_FIGS
    % Create the directory (if it doesn't already exist)
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end

    % Save the figure in the results folder
    saveas(gcf, fullfile(rf_fldr, fil_nm));
    cd ..
    cd ..
    cd Matlab_codes

    %% Generating and plotting the comparison in frequency domain
    figure
    for i_s = 1:n_c

        for i_col = 1:size(sig_mat{i_s}, 2)
            sig_vect = sig_mat{i_s}(:, i_col);
            %%
            Fs=1./dt; % sampling rate
            n_fft=length(t_vect); %number of samples
            [sig_freq,freq_ss,~]= fns_data_process.fun_fftandifft(t_vect,Fs,sig_vect);
            %%
            t_seissol=data_seissol.t;
            data_vect=data_seissol{:,i_s+1};
            Fs_seissol=1./(t_seissol(3)-t_seissol(2)); % sampling rate
            n_fft=length(t_seissol); %number of samples
            [data_seisol_freq,freq2,~]= fns_data_process.fun_fftandifft(t_seissol,Fs_seissol,data_vect);
            %%
            % Plot the time-domain signal
            %         figure
            l=l_vect(i_col);
            b=b_vect(i_col);
            txt = ['SeisSol';'LPM-FEM'];
            %             txt = 'LPM-FEM';
            txt_2 = ['Floor:',num2str(i_str),...
                ',~Component:',num2str(i_s)];
            subplot(n_c, size(sig_mat{i_s}, 2),...
                (i_s-1)*size(sig_mat{i_s}, 2) + i_col);
            plot(freq2, abs(data_seisol_freq), 'LineWidth', 1.5)
            hold on
            plot(freq_ss, abs(sig_freq),'LineStyle',':','LineWidth', 1.5);
            xlabel('frequency (hz)','Interpreter','latex','FontSize',12);
            ylabel('Velocity (m/s)','Interpreter','latex','FontSize',12);
            xlim([0, 40])
            title(txt_2,'Interpreter','latex','FontSize',11);
            legend(txt);
            legend('Box','off','Interpreter','latex','FontSize',11)
            set(gcf,'Units','inches', 'Position', [18 3 4.5 6],...
                'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125]);
        end
    end
    fil_nm = ['Fdomain_vel_t_','_Floor_', num2str(i_str),...
        '_ftyp_', ftyp,'_Vs_', num2str(V_s),...
        '_Lf_', num2str(L_f), '_Bf_', num2str(B_f), '.emf'];
    cd SAVE_FIGS
    % Create the directory (if it doesn't already exist)
    if ~exist(rf_fldr, 'dir')
        mkdir(rf_fldr);
    end

    % Save the figure in the results folder
    saveas(gcf, fullfile(rf_fldr, fil_nm));
    cd ..
    cd ..
    cd Matlab_codes
end
