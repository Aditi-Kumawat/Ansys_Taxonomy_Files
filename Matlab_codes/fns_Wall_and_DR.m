classdef fns_Wall_and_DR
    methods (Static)
        %%
        function [f_vect,TFamp_mat,TFcmplx_mat]=...
                get_TFWall(n_str, n_rx, n_ry,l, b,...
                ftyp, V_s, L_f, B_f,array_para,bf_nm,i_str,cmpt,n_c,r_fldr,cols)

            for i_para = 1:length(array_para)
                para=array_para(i_para);
                % Get the folder name for the specific configuration
                fldr = fns_Wall_and_DR.get_fldrnmWall(n_str,...
                    n_rx, n_ry, l, b, ftyp, V_s, L_f, B_f, num2str(para));
                fl_nm1 = arrayfun(@(x) sprintf(bf_nm, x{1}, i_str,...
                    l, b),cmpt, 'UniformOutput', false);
                cd ..
                cd Results_Ansys
                fil_pth = fullfile(r_fldr, fldr, fl_nm1);
                U_all = cellfun(@(x) readtable(x), fil_pth,...
                    'UniformOutput', false);
                cd ..
                cd Matlab_codes

                for i_cmp = 1:n_c
                    U = U_all{i_cmp};
                    U.Properties.VariableNames = cols;
                    if i_para == 1
                        f_vect = U.Freq;
                        TFamp_mat{i_cmp} = U.AMPL;
                        TFr_mat{i_cmp} = U.REAL;
                        TFim_mat{i_cmp} = U.IMAG;
                        TFcmplx_mat{i_cmp} = U.REAL+1i.*U.IMAG;
                    else
                        TFamp_mat{i_cmp}(:, i_para) = U.AMPL;
                        TFr_mat{i_cmp}(:, i_para)  = U.REAL;
                        TFim_mat{i_cmp}(:, i_para)  = U.IMAG;
                        TFcmplx_mat{i_cmp}(:, i_para)  = ...
                            U.REAL+1i.*U.IMAG;
                    end
                end
            end
        end
        %%
        function [f_vect,TFamp_mat,TFcmplx_mat]=...
                get_TF_DR(n_str, n_rx, n_ry,l, b,...
                ftyp, V_s, L_f, B_f,array_para,bf_nm,i_str,cmpt,n_c,r_fldr,cols)

            for i_para = 1:length(array_para)
                para=array_para(i_para);
                % Get the folder name for the specific configuration
                fldr = fns_Wall_and_DR.get_fldrnmDR(n_str,...
                    n_rx, n_ry, l, b, ftyp, V_s, L_f, B_f, num2str(para));
                fl_nm1 = arrayfun(@(x) sprintf(bf_nm, x{1}, i_str,...
                    l, b),cmpt, 'UniformOutput', false);
                cd ..
                cd Results_Ansys
                fil_pth = fullfile(r_fldr, fldr, fl_nm1);
                U_all = cellfun(@(x) readtable(x), fil_pth,...
                    'UniformOutput', false);
                cd ..
                cd Matlab_codes

                for i_cmp = 1:n_c
                    U = U_all{i_cmp};
                    U.Properties.VariableNames = cols;
                    if i_para == 1
                        f_vect = U.Freq;
                        TFamp_mat{i_cmp} = U.AMPL;
                        TFr_mat{i_cmp} = U.REAL;
                        TFim_mat{i_cmp} = U.IMAG;
                        TFcmplx_mat{i_cmp} = U.REAL+1i.*U.IMAG;
                    else
                        TFamp_mat{i_cmp}(:, i_para) = U.AMPL;
                        TFr_mat{i_cmp}(:, i_para)  = U.REAL;
                        TFim_mat{i_cmp}(:, i_para)  = U.IMAG;
                        TFcmplx_mat{i_cmp}(:, i_para)  = ...
                            U.REAL+1i.*U.IMAG;
                    end
                end
            end
        end
        %%
        function fldr_nm = get_fldrnmWall(n_str, n_rx, n_ry, l, b,...
                ftyp,V_s, L_f, B_f,para)
            fldr_nm = ['n_storeys_',num2str(n_str),'_n_rooms_X_',...
                num2str(n_rx),'_n_rooms_Y_',num2str(n_ry),...
                '_l',num2str(l),'_by_b',num2str(b),'_ftyp_',ftyp,...
                '_Vs_',num2str(V_s),...
                '_Lf_',num2str(L_f),'_Bf_',num2str(B_f),...
                '_CONFIG_',para];
        end
        %%
        function fldr_nm = get_fldrnmDR(n_str, n_rx, n_ry, l, b,...
                ftyp,V_s, L_f, B_f,para)
            fldr_nm = ['n_storeys_',num2str(n_str),'_n_rooms_X_',...
                num2str(n_rx),'_n_rooms_Y_',num2str(n_ry),...
                '_l',num2str(l),'_by_b',num2str(b),'_ftyp_',ftyp,...
                '_Vs_',num2str(V_s),...
                '_Lf_',num2str(L_f),'_Bf_',num2str(B_f),...
                '_DR_',para];
        end
        %%
        function [f_vect, TF_amp_mat, TF_cpmlx_mat] = get_TF(rf_fldr, n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f, wall_config, dampg_vect, i_str, n_c)
            bf_nm = 'Disp_Center_%s_%d_l%d_b%d';
            cols1 = {'Freq', 'AMPL','PHASE','REAL','IMAG'};
            cmpt = {'X', 'Y', 'Z'};
            % Initialize bld_cases to an empty array for the cases where it's not used
            if strcmp(rf_fldr, 'Bld_with_Walls')
                [f_vect, TF_amp_mat, TF_cpmlx_mat] = fns_Wall_and_DR.get_TFWall(n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f, wall_config, bf_nm, i_str, cmpt, n_c, rf_fldr, cols1);
            elseif strcmp(rf_fldr, 'Vary_DampRatio')
                [f_vect, TF_amp_mat, TF_cpmlx_mat] = fns_Wall_and_DR.get_TF_DR(n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f, dampg_vect, bf_nm, i_str, cmpt, n_c, rf_fldr, cols1);
            else
                [f_vect, TF_amp_mat, TF_cpmlx_mat, lb_cmbs] = fns_scatter.get_TF_scatter(n_str, n_rx, n_ry, l_vect, b_vect, ftyp, V_s, L_f, B_f, bf_nm, i_str, cmpt, n_c, rf_fldr, cols1);

            end
        end
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

    end
end