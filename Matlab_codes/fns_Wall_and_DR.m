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
    end
end