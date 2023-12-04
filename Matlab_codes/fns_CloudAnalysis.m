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

    end
end