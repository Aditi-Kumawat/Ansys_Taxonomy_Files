classdef fns_CloudAnalysis
    methods (Static)
        %%
        function save_DINvals(data_set,date_evnt,time_evnt,...
                rf_fldr,bld_soil_fndn,f_x_max,f_y_max,f_z_max,...
                max_Vxmat,max_Vymat,max_Vzmat,t_x_max,t_y_max,t_z_max,...
                max_Vx_KB_f_mat,max_Vy_KB_f_mat,max_Vz_KB_f_mat)
            fldr_evnt=['DINvals_', date_evnt, '_', time_evnt];
            cd SAVE_DATA;
            if ~exist(data_set, 'dir')
                mkdir(data_set);
            end
            cd(data_set)
            if ~exist(rf_fldr, 'dir')
                mkdir(rf_fldr);
            end
            cd(rf_fldr)
            if ~exist(bld_soil_fndn, 'dir')
                mkdir(bld_soil_fndn);
            end
            cd(bld_soil_fndn);
            if ~exist(fldr_evnt, 'dir')
                mkdir(fldr_evnt);
            end
            cd(fldr_evnt);
            % Specify the file name for the .mat file
            filename1 = 'DIN4150_3_Vmax.mat';
            filename2 = 'DIN4150_2_KBval.mat';
            % Save the variables to the .mat file
            save(filename1, 'f_x_max', 'f_y_max', 'f_z_max', ...
                'max_Vxmat', 'max_Vymat', 'max_Vzmat');
            save(filename2, 't_x_max', 't_y_max', 't_z_max', ...
                'max_Vx_KB_f_mat', 'max_Vy_KB_f_mat', 'max_Vz_KB_f_mat');
            cd(fullfile('..', '..', '..', '..', '..'));
        end
        %%
        function [f_x_max, f_y_max, f_z_max, max_Vxmat, max_Vymat, max_Vzmat, ...
                t_x_max, t_y_max, t_z_max, max_Vx_KB_f_mat, max_Vy_KB_f_mat, max_Vz_KB_f_mat] = ...
                import_DINvals(data_set, date_evnt, time_evnt, rf_fldr, bld_soil_fndn)

            % Construct the directory path
            baseDir = 'SAVE_DATA'; % Adjust this if your base directory is different
            fldr_evnt = ['DINvals_', date_evnt, '_', time_evnt];
            dirPath = fullfile(baseDir, data_set, rf_fldr, bld_soil_fndn, fldr_evnt);

            % Check if the directory exists
            if ~exist(dirPath, 'dir')
                error('Directory %s does not exist.', dirPath);
            end

            % Define file names
            filename1 = fullfile(dirPath, 'DIN4150_3_Vmax.mat');
            filename2 = fullfile(dirPath, 'DIN4150_2_KBval.mat');

            % Load data from the files
            if exist(filename1, 'file')
                data1 = load(filename1);
                f_x_max = data1.f_x_max;
                f_y_max = data1.f_y_max;
                f_z_max = data1.f_z_max;
                max_Vxmat = data1.max_Vxmat;
                max_Vymat = data1.max_Vymat;
                max_Vzmat = data1.max_Vzmat;
            else
                error('File %s does not exist.', filename1);
            end

            if exist(filename2, 'file')
                data2 = load(filename2);
                t_x_max = data2.t_x_max;
                t_y_max = data2.t_y_max;
                t_z_max = data2.t_z_max;
                max_Vx_KB_f_mat = data2.max_Vx_KB_f_mat;
                max_Vy_KB_f_mat = data2.max_Vy_KB_f_mat;
                max_Vz_KB_f_mat = data2.max_Vz_KB_f_mat;
            else
                error('File %s does not exist.', filename2);
            end
        end
                %%
        function [max_Vxmat, max_Vymat, max_Vzmat,...
                max_Vx_KB_f_mat, max_Vy_KB_f_mat, max_Vz_KB_f_mat] = ...
                import_DINvalsSingleBld(data_set, date_evnt, time_evnt, rf_fldr, bld_soil_fndn)

            % Construct the directory path
            baseDir = 'SAVE_DATA'; % Adjust this if your base directory is different
            fldr_evnt = ['DINvals_', date_evnt, '_', time_evnt];
            dirPath = fullfile(baseDir, data_set, rf_fldr, bld_soil_fndn, fldr_evnt);

            % Check if the directory exists
            if ~exist(dirPath, 'dir')
                error('Directory %s does not exist.', dirPath);
            end

            % Define file names
            filename1 = fullfile(dirPath, 'singlebldDIN4150_3_Vmax.mat');
            filename2 = fullfile(dirPath, 'singlebldDIN4150_2_KBval.mat');

            % Load data from the files
            if exist(filename1, 'file')
                data1 = load(filename1);
                max_Vxmat = data1.max_Vxmat;
                max_Vymat = data1.max_Vymat;
                max_Vzmat = data1.max_Vzmat;
            else
                error('File %s does not exist.', filename1);
            end

            if exist(filename2, 'file')
                data2 = load(filename2);
                max_Vx_KB_f_mat = data2.max_Vx_KB_f_mat;
                max_Vy_KB_f_mat = data2.max_Vy_KB_f_mat;
                max_Vz_KB_f_mat = data2.max_Vz_KB_f_mat;
            else
                error('File %s does not exist.', filename2);
            end
        end
        %%
        function save_freq_forVmax(ff_Vamp_mat,f_inpt_V,ff_fldr)
            maxVx=max(ff_Vamp_mat{1});
            maxVy=max(ff_Vamp_mat{2});
            maxVz=max(ff_Vamp_mat{3});
            [maxrowidx_vx, ~] = find(bsxfun(@eq, ff_Vamp_mat{1}, maxVx));
            [maxrowidx_vy, ~] = find(bsxfun(@eq, ff_Vamp_mat{2}, maxVy));
            [maxrowidx_vz, ~] = find(bsxfun(@eq, ff_Vamp_mat{3}, maxVz));
            fx=f_inpt_V{1};
            fy=f_inpt_V{2};
            fz=f_inpt_V{3};
            f_vx_max=fx(maxrowidx_vx);
            f_vy_max=fy(maxrowidx_vy);
            f_vz_max=fz(maxrowidx_vz);
            filename1 = 'freq_Hz_max_vamp.mat';
            cd ..
            fullFilePath = fullfile(ff_fldr, filename1);
            save(fullFilePath, 'f_vx_max', 'f_vy_max', 'f_vz_max');
            cd Matlab_codes
        end
        %%
        function save_PGV(ff_Vt,ff_fldr)
            PGVx=max(abs(ff_Vt{1}));
            PGVy=max(abs(ff_Vt{2}));
            PGVz=max(abs(ff_Vt{3}));
            filename2 = 'PGV_xyz.mat';
            cd ..
            fullFilePath = fullfile(ff_fldr, filename2);
            save(fullFilePath, 'PGVx', 'PGVy', 'PGVz');
            cd Matlab_codes
        end
        %%
        function [PGVx, PGVy, PGVz,fx,fy,fz] = import_PGV(ff_fldr)
            % Import PGV data from a specified folder
            %
            % Args:
            %   ff_fldr (string): Path to the folder containing the PGV data file
            %
            % Returns:
            %   PGVx, PGVy, PGVz (arrays): Peak Ground Velocity in x, y, z directions

            % Construct the full path to the file
            filename1 = 'freq_Hz_max_vamp.mat';
            filename2 = 'PGV_xyz.mat';
            fullFilePath1 = fullfile(ff_fldr, filename1);
            fullFilePath2 = fullfile(ff_fldr, filename2);
            cd ..
            % Load the data from the file
            data1 = load(fullFilePath1);
            data2 = load(fullFilePath2);

            % Extract PGV values
            PGVx = data2.PGVx;
            PGVy = data2.PGVy;
            PGVz = data2.PGVz;
            % Extract fmax values
            fx = data1.f_vx_max;
            fy = data1.f_vy_max;
            fz = data1.f_vz_max;
            cd Matlab_codes
        end
        %% plot results
        function plot_CA(event_vect,all_PGV,all_max_v_or_KB,event_labels,ylb,xlb,titleText,filename)
            evntlbl_length = jet(length(event_vect));
            figure;
            fnt=16;
            s = scatter(all_PGV, all_max_v_or_KB, [], event_labels, 'filled');
            colormap(evntlbl_length);
            s.MarkerEdgeColor = 'k';
            s.LineWidth = 0.1;
            s.MarkerFaceAlpha = 1;
            % Calculate the maximum value of all_PGV and add some space
            max_PGV = max(all_PGV);
            space_after_max = max_PGV * 0.1; % for example, 10% more than the max value
            xlim([min(all_PGV) max_PGV + space_after_max]); % Set new x-axis limits
            grid on;
            xlabel(xlb, 'FontSize', fnt,'Interpreter', 'latex');
            ylabel(ylb, 'FontSize', fnt,'Interpreter', 'latex');
            set(gca, 'XScale', 'log');
            set(gca, 'FontSize', fnt, 'Box', 'on', 'LineWidth', 0.5,...
                'TickLabelInterpreter', 'latex', 'TickLength',[0.01,0.01]);
            set(gcf, 'Units', 'inches', 'Position', [18 3 6 4],...
                'PaperUnits', 'Inches', 'PaperSize', [6 4]);
            % Create a custom legend
            hold on;
            legend_entries = gobjects(length(event_vect), 1); % Initialize an array of graphic objects
            for ie = 1:length(event_vect)
                % Create scatter plots for legend entries with the same marker properties
                legend_entries(ie) = scatter(NaN, NaN, [], evntlbl_length(ie, :), 'filled', 'DisplayName', event_vect{ie});
                legend_entries(ie).MarkerEdgeColor = 'k';
                legend_entries(ie).LineWidth = 0.1;
                legend_entries(ie).MarkerFaceAlpha = 1;
            end

            % Set legend to two columns
            legend(legend_entries, 'NumColumns', 3, 'Location', 'northwest', 'FontSize', fnt,...
                'Interpreter', 'latex','Box','off');
            title(titleText, 'Interpreter', 'latex', 'FontSize', fnt);
            hold off;
            fns_FragilityFns.save_figure_to_pdf(filename);
        end
        %% plot results
        function plot_CAfull(event_vect,all_PGV,all_max_v_or_KB,event_labels,ylb,xlb,filename,ylimits)
            evntlbl_length = jet(length(event_vect));

            colors = evntlbl_length(event_labels, :);
            figure;
            fnt=12;
            for i = 1:size(all_max_v_or_KB, 2)
                s = scatter(all_PGV, all_max_v_or_KB(:,i), 50, colors, 'filled');
                colormap(evntlbl_length);
                s.MarkerEdgeColor = 'k';
                s.LineWidth = 0.1;
                % transperency commands-don't work with legend
                % s.MarkerFaceAlpha = 0.6;
                % s.MarkerEdgeAlpha = 0.6;
                hold on;  % Hold on to plot all columns on the same figure
            end
            % Calculate the maximum value of all_PGV and add some space
            max_PGV = max(all_PGV);
            space_after_max = max_PGV * 0.5; % for example, 10% more than the max value
            xlim([min(all_PGV) max_PGV + space_after_max]); % Set new x-axis limits
            ylim(ylimits);
            grid on;
            xlabel(xlb, 'FontSize', fnt,'Interpreter', 'latex');
            ylabel(ylb, 'FontSize', fnt,'Interpreter', 'latex');
            set(gca, 'XScale', 'log');
            set(gca, 'FontSize', fnt, 'Box', 'on', 'LineWidth', 0.5,...
                'TickLabelInterpreter', 'latex', 'TickLength',[0.01,0.01]);
            set(gcf, 'Units', 'inches', 'Position', [18 3 5 3.5],...
                'PaperUnits', 'Inches', 'PaperSize', [5 3.5]);
            % Create a custom legend
            hold on;
            legend_entries = gobjects(length(event_vect), 1);  % Initialize an array of graphic objects
            for ie = 1:length(event_vect)
                % Create dummy scatter plot handles for legend entries
                legend_entries(ie) = plot(NaN, NaN, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', evntlbl_length(ie, :), 'DisplayName', event_vect{ie});
                % set(legend_entries(ie), 'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6);
            end

            % Set legend with the custom plot handles
            legend(legend_entries, 'NumColumns', 2, 'Location', 'northwest', 'FontSize', fnt, 'Interpreter', 'latex', 'Box', 'off');

            % title(titleText, 'Interpreter', 'latex', 'FontSize', fnt);
            hold off;
            fns_FragilityFns.save_figure_to_pdf(filename);
        end
    end
end