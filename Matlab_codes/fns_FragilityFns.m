classdef fns_FragilityFns
    methods (Static)
        function [sorted_PGV,sorted_fnvals,fitted_fn]=...
                get_curve_polyft(evry_maxval,all_PGV,lim,ndeg)
            fraction_lowlimKB_x = zeros(size(evry_maxval, 1), 1);
            % Calculate the fraction for each PGV value
            for i = 1:size(evry_maxval, 1)
                fraction_lowlimKB_x(i) =...
                    sum(evry_maxval(i, :) >= lim)/size(evry_maxval, 2);
            end
            [sorted_PGV, idx_x] = sort(all_PGV);
            sorted_fnvals = fraction_lowlimKB_x(idx_x);
            p_x = polyfit(sorted_PGV, sorted_fnvals, ndeg);
            fitted_fn = polyval(p_x, sorted_PGV);
        end
        %%
        function [sorted_PGV, sorted_fnvals, fitted_fn] =...
                get_curve_splineft(evry_maxval, all_PGV, lim)
            fraction_lowlimKB_x = zeros(size(evry_maxval, 1), 1);
            % Calculate the fraction for each PGV value
            for i = 1:size(evry_maxval, 1)
                fraction_lowlimKB_x(i) = ...
                    sum(evry_maxval(i, :) > lim)/size(evry_maxval, 2);
            end
            [sorted_PGV, idx_x] = sort(all_PGV);
            sorted_fnvals = fraction_lowlimKB_x(idx_x);
            % Spline fitting
            spline_fit = spline(sorted_PGV, sorted_fnvals);
            % Evaluate the spline at the sorted PGV values
            fitted_fn = ppval(spline_fit, sorted_PGV);
        end
        %%
        function [sorted_PGV, sorted_fnvals, fitted_fn] =...
                get_curve_splineft2(evry_maxval, all_PGV, lim, smoothingParam)
            fraction_lowlimKB_x = zeros(size(evry_maxval, 1), 1);
            % Calculate the fraction for each PGV value
            for i = 1:size(evry_maxval, 1)
                fraction_lowlimKB_x(i) = ...
                    sum(evry_maxval(i, :) > lim)/size(evry_maxval, 2);
            end
            [sorted_PGV, idx_x] = sort(all_PGV);
            sorted_fnvals = fraction_lowlimKB_x(idx_x);
            % Spline fitting with smoothing parameter
            spline_fit = csaps(sorted_PGV, sorted_fnvals, smoothingParam);
            % Evaluate the spline at the sorted PGV values
            fitted_fn = ppval(spline_fit, sorted_PGV);
            % Post-process to ensure no dips below the peak
            % (assuming the peak is at 1)
            fitted_fn(fitted_fn > 1) = 1;
        end

        %%
        function plot_direction(direction, limKB_vect, evry_max_KB,...
                all_PGV, ndeg, evry_max_V, lim_v, lcol, legend_vect)
            % Decide which axis label to use based on the direction
            if strcmp(direction, 'x')
                axisLabel = 'log(PGV($x$)), m/s';
                titleText = 'CDF for Building Responses ($x$-dir)';
                filename = 'CDFPlot_x.pdf';
            elseif strcmp(direction, 'y')
                axisLabel = 'log(PGV($y$)), m/s';
                titleText = 'CDF for Building Responses ($y$-dir)';
                filename = 'CDFPlot_y.pdf';
            else
                axisLabel = 'log(PGV($z$)), m/s';
                titleText = 'CDF for Building Responses ($z$-dir)';
                filename = 'CDFPlot_z.pdf';
            end
            fnt=14;
            figure
            for i = 1:length(limKB_vect)
                limKB = limKB_vect(i);
                [sorted_PGV, sorted_P_KB, fitted_P_KB] =...
                    fns_FragilityFns.get_curve_polyft(evry_max_KB,...
                    all_PGV, limKB, ndeg);
                % Ensure fitted values are within the [0,1] range
                fitted_P_KB(fitted_P_KB < 0) = 0;
                fitted_P_KB(fitted_P_KB > 1) = 1;
                % Find the maximum value in the fitted curve
                max_fitted_value = max(fitted_P_KB);
                % Find the index of the first occurrence of the maximum value
                idx_max_value = find(fitted_P_KB == max_fitted_value, 1, 'first');
                % Ensure that all subsequent values are set to the maximum value
                fitted_P_KB(idx_max_value:end) = max_fitted_value;

                % scatter(sorted_PGV, sorted_P_KB, 'linewidth', 1.2,...
                %     'MarkerEdgeColor', lcol{i})
                hold on;
                plot(sorted_PGV, fitted_P_KB, 'linewidth', 2,...
                    'Color', lcol{i})
            end

            [~, sorted_P_V, fitted_P_V] =...
                fns_FragilityFns.get_curve_polyft(evry_max_V,...
                all_PGV,lim_v, ndeg);
            % scatter(sorted_PGV, sorted_P_V, 'linewidth', 1.2,...
            %     'MarkerEdgeColor', lcol{end})
            plot(sorted_PGV, fitted_P_V,':', 'linewidth', 1.8,...
                'Color', lcol{end})
            ylim([0 1])
            xlim([min(sorted_PGV) max(sorted_PGV)]);
            set(gca, 'XScale', 'log');
            xlabel(axisLabel, 'Interpreter', 'latex', 'FontSize', fnt);
            ylabel('Prob($v_{i} \geq lim |$PGV)', 'Interpreter',...
                'latex', 'FontSize', fnt);
            legend(legend_vect, 'Interpreter', 'latex', 'Location',...
                'northwest', 'Box', 'off');
            % title(titleText, 'Interpreter', 'latex', 'FontSize', fnt);
            grid on;
            set(gca, 'FontSize', fnt, 'Box', 'on', 'LineWidth', 0.5,...
                'TickLabelInterpreter', 'latex', 'TickLength',...
                [0.01, 0.01]);

            set(gcf, 'Units', 'inches', 'Position', [18 3 5.5 3.5],...
                'PaperUnits', 'Inches', 'PaperSize', [5.5 3.5]);
            fns_FragilityFns.save_figure_to_pdf(filename);
        end

        function save_figure_to_pdf(filename)
            cd SAVE_FIGS
            cd Fragility_Fn
            saveas(gcf, filename);
            cd ..
            cd ..
        end

    end
end