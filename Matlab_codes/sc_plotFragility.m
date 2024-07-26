% This MATLAB script loads fragility data from a .mat file and plots fragility functions.
% The fragility functions are estimated using Maximum Likelihood Estimation (MLE) method.
%% import sorted_PGV, KB_exceed_mat, v_exceed_mat, nbld
% saved using sc_FragilityFn.m
close all; clear; clc
cd SAVE_DATA
cd Fragility_Fn
load fragilitydata.mat
cd ..
cd ..
%% assign variables for the code by Baker 2015
IM = sorted_PGV*1e3;
num_gms = nbld * ones(1, length(sorted_PGV(1, :)));
x_vals = (1e-6:1e-6:1e0)*1e3; % IM levels to plot fragility function at
%% plotting
axisLabel = {'log(PGV($x$)), mm/s', 'log(PGV($y$)), mm/s', 'log(PGV($z$)), mm/s'};
filename = 'CDFPlot_x.pdf';

plotFragility(KB_exceed_mat, v_exceed_mat(1, :), IM(1, :), num_gms, x_vals, axisLabel{1}, filename)
% plotFragility(KB_exceed_mat, v_exceed_mat(2, :), IM(2, :), num_gms, x_vals, axisLabel{2})

%% ------------------------------------------------------------------------------------%
function plotFragility(KB_exceed_mat, v_exceed_vect, IM, num_gms, x_vals, axisLabel, filename)
% the function computes the fragility functions and plots them with PGV (IM),
% the computation differs from the usual fragility computations as here we
% have the num_gms that are equivalent to the number of test cases for the
% building model and not equal to the number of ground motions
% corresponding to the same IM. 
    
    ha_cl = @colors;
    lcol = {ha_cl('cinnamon'), ha_cl('crimson'), ...
            ha_cl('gray'), ha_cl('denim'), ha_cl('black')};
    fnt = 11;
    figure
    legend_handles = []; % Array to store handles for the legend
    %% fragility curve for velocity value
    sorted_collapse = sort(v_exceed_vect);
    num_collapse = sorted_collapse;
    % estimate fragility function using MLE method (equation 11)
    [theta_hat_mle, beta_hat_mle] = fns_FragilityFns.fn_mle_pc(IM, num_gms, num_collapse);
    % compute fragility function using equation 1 and estimated parameters
    p_collapse_mle = normcdf((log(x_vals / theta_hat_mle)) / beta_hat_mle);
    h2 = plot(IM, num_collapse ./ num_gms, '^', 'linewidth', 1, 'Color', lcol{end});
    legend_handles = [legend_handles, h2]; % Store handle for legend
    hold on
    %% fragility curve for KB values
    plot(x_vals, p_collapse_mle, '--', 'linewidth', 1.5, 'Color', lcol{end})
    for i = 1:size(KB_exceed_mat, 1)
        sorted_collapse = sort(KB_exceed_mat(i, :));
        num_collapse = sorted_collapse;
        % estimate fragility function using MLE method (equation 11)
        [theta_hat_mle, beta_hat_mle] = fns_FragilityFns.fn_mle_pc(IM, num_gms, num_collapse);
        % compute fragility function using equation 1 and estimated parameters
        p_collapse_mle = normcdf((log(x_vals / theta_hat_mle)) / beta_hat_mle);
        h1 = plot(IM, num_collapse ./ num_gms, '^', 'linewidth', 1, 'Color', lcol{i});
        legend_handles = [legend_handles, h1]; % Store handle for legend
        hold on
        plot(x_vals, p_collapse_mle, '-', 'linewidth', 1.5, 'Color', lcol{i})
        hold on
    end
    %%
    xlim([min(x_vals) max(x_vals)]);
    ylim([-0.1 1.1]);
    set(gca, 'XScale', 'log');

    xlabel(axisLabel, 'Interpreter', 'latex', 'FontSize', fnt);
    ylabel('P($v \geq lim |$PGV($x$))', 'Interpreter', ...
           'latex', 'FontSize', fnt);
    legend_vect = {'$\max[v_x(t)] \geq v_x(lim)$'...
                   '$\max[KB_{F}] \geq A_u(d)$', ...
                   '$\max[KB_{F}] \geq A_0(d)$', ...
                   '$\max[KB_{F}] \geq A_u(n)$', ...
                   '$\max[KB_{F}] \geq A_0(n)$'};
    % legend_vect = {'Frac($\max[KB_{F}] \geq A_u(d)$)', ...
    %                '$\max[KB_{F}] \geq A_o(d)$', ...
    %                '$\max[KB_{F}] \geq A_u(n)$', ...
    %                '$\max[KB_{F}] \geq A_o(n)$', ...
    %                '$\max[v_x(t)] \geq v_x(lim)$'};
    legend(legend_handles, legend_vect, 'FontSize', fnt, 'Interpreter',...
        'latex', 'Location', 'east', 'Box', 'off');
    % title(titleText, 'Interpreter', 'latex', 'FontSize', fnt);
    grid on;
    set(gca, 'FontSize', fnt, 'Box', 'on', 'LineWidth', 0.5, ...
             'TickLabelInterpreter', 'latex', 'TickLength', ...
             [0.01, 0.01]);

    set(gcf, 'Units', 'inches', 'Position', [18 3 9 2.5], ...
             'PaperUnits', 'Inches', 'PaperSize', [9 2.5]);
    %% save figure
    fns_FragilityFns.save_figure_to_pdf(filename);
end
