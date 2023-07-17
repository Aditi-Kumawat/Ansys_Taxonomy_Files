
%% Initialization
clear;clc;close all;

% Define the number of storeys, rooms in x- y-direction
n_str = 1;
n_rx = 1;
n_ry = 1;

% Define the length, width, and height of the building
l_vect=[2 3 4 5 6 7 8];
b_vect=[2 3 4 5 6 7 8];
% l_vect=[3 5 7];
% b_vect=[3 5 7];
h = 3;

% Define the type of foundation as either 'PLATE' or 'FOOTING'
ftyp = 'FOOTING';

% Define the velocity of the excitation
V_s = 450;

% Define the size of the elements
n_esize = 0.25;

% Calculate the length and width of the footing based on the
% foundation type
if strcmp(ftyp,'PLATE')
    B_f = n_esize/2;
    L_f = n_esize/2;
else
    B_f = 0.75;
    L_f = 0.75;
end
%% Importing Transfer Function
% Define the name of the folder where the results are stored
rf_fldr = 'UnitBld_GeomVary';
bf_nm = 'Disp_Center_%s_%d_l%d_b%d';

% Define the column names in the results table
cols = {'Freq', 'AMPL','PHASE','REAL','IMAG'};

% Define the components to extract from the results
cmpt = {'X', 'Y', 'Z'};
%%
n_c = length(cmpt);
Vabs_zCell = cell(1, n_str+1);
fref_cmplx_mat_1=cell(1, n_c);
DISP_cmplx_mat=cell(1, n_c);
VEL_abs_mat=cell(1, n_c);
fref_vel_amp_mat_1=cell(1, n_c);

for i_str = 0:n_str
    [f_vect,TF_amp_mat,TF_cpmlx_mat,lb_combs]=...
        fns_scatter.get_TF_scatter(n_str, n_rx, n_ry,...
        l_vect, b_vect, ftyp, V_s, L_f, B_f,...
        bf_nm,i_str,cmpt,n_c,rf_fldr,cols);

    %%
    TFabs_zcell{i_str+1} = TF_amp_mat{3};
    TFabs_xcell{i_str+1} = TF_amp_mat{1};
end
%%
plt_cmp='Z';
y_lim=50;
TF_abs_Zmat=TFabs_zcell{n_str+1};
n_blds = size(TF_abs_Zmat, 2);
[max_vals, max_indices] = max(TF_abs_Zmat);
f_max_vals = f_vect(max_indices);
max_vals_all = [];
f_max_vals_all = [];
max_vals_all = [max_vals_all, max_vals];
f_max_vals_all = [f_max_vals_all, f_max_vals];

% Convert lb_combs to string format
lb_combs_str = cellstr(num2str(lb_combs));
% Convert lb_combs to string format with the desired label format
lb_combs_str1 = cell(size(lb_combs, 1), 1);
for i = 1:size(lb_combs, 1)
    lb_combs_str1{i} = sprintf('%dmx%dm', lb_combs(i, 1), lb_combs(i, 2));
end
% Get unique values from lb_combs_str
unique_lb_combs_str = unique(lb_combs_str);
unique_lb_combs_str1 = unique(lb_combs_str1);

% Determine the corresponding color indices for each unique lb_combs_str value
color_indices = zeros(size(lb_combs_str));
for i = 1:numel(unique_lb_combs_str)
    idx = strcmp(lb_combs_str, unique_lb_combs_str{i});
    color_indices(idx) = i;
end
% create horizontal line for the scatter plot
y = ones(size(max_vals_all));

% create colormap based on maximum values
c = max_vals_all;
cmap = jet(length(c));
% Create a discrete colormap based on unique_lb_combs_str
discrete_cmap = colormap(jet(numel(unique_lb_combs_str)));

% Assign colors based on color_indices
c = color_indices;

% Make scatter plot with discrete colormap
s = scatter(f_max_vals_all, max_vals_all, 100, c, 'filled');
colormap(discrete_cmap);
s.MarkerEdgeColor = 'k';
s.LineWidth = 1;
s.MarkerFaceAlpha = 0.8;
% Add colorbar with unique_lb_combs_str as tick labels
cbar = colorbar('Ticks', 1:numel(unique_lb_combs_str1), 'TickLabels', unique_lb_combs_str1);

%             cbar.Label.String = 'lb_combs';
cbar.Label.FontSize = 12;
cbar.TickLabelInterpreter = 'latex';

% adjust plot properties
xlim([min(f_vect), y_lim]);
ylim([2, 12]);
xlabel('Frequency (Hz)', 'FontSize', 14,...
    'Interpreter', 'latex');
ylabel('Maximum values', 'FontSize', 14,...
    'Interpreter', 'latex');
%             title('Scatter plot of Maximum values vs Frequency',...
%                 'FontSize', 14, 'FontWeight', 'bold',...
%                 'Interpreter', 'latex');
set(gca, 'FontSize', 14, 'Box', 'on', 'LineWidth', 1,...
    'TickLabelInterpreter', 'latex', 'TickLength',[0.01,0.01]);
set(gcf, 'Units', 'inches', 'Position', [18 3 10 6],...
    'PaperUnits', 'Inches', 'PaperSize', [7.25, 9.125]);
cbar.Label.FontSize = 12;

filnm = ['TF_Scatter_Plot_',plt_cmp, '.emf'];

cd SAVE_FIGS
if ~exist(rf_fldr, 'dir')
    mkdir(rf_fldr);
end
saveas(gcf, fullfile(rf_fldr, filnm));
cd ..
cd ..
cd Matlab_codes