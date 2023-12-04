%% fns_Inpt_Para
classdef fns_Inpt_BldPara
    methods (Static)
        %%
        function selectedOption = selectRfFldr()
            options = {'MultiUnitBld_GeomVary_3lby2', 'MultiUnitBld_GeomVary', 'Bld_with_Walls', 'Vary_DampRatio'};
            [choice, isOk] = listdlg('PromptString', 'Select an option:', ...
                'SelectionMode', 'single', ...
                'ListString', options,'ListSize',[350,150]);

            if isOk
                selectedOption = options{choice};
            else
                selectedOption = 'None';
                disp('No valid selection made.');
            end
        end

        %%
        function selectedOption = select_bldsoilpara()
            options = {'nstr3_Plate450', 'nstr3_Plate100', 'nstr1_Footing450', 'nstr1_Footing100', 'nstr2_Plate450', 'nstr2_Plate100'};
            [choice, isOk] = listdlg('PromptString', 'Select an option:', ...
                'SelectionMode', 'single', ...
                'ListString', options,'ListSize',[150,150]);

            if isOk
                selectedOption = options{choice};
            else
                selectedOption = 'None';
                disp('No valid selection made.');
            end
        end

        %%
        function [l_vect,b_vect,h,wall_config,dampg_vect,bld_cases]=get_lbh_bldcases_for_rf_fldr(rf_fldr)
            if strcmp(rf_fldr, 'MultiUnitBld_GeomVary_3lby2') || strcmp(rf_fldr, 'MultiUnitBld_GeomVary')
                l_vect=[2 3 4 5 6 7 8];
                b_vect=[2 3 4 5 6 7 8];
                h = 3;
                wall_config=[];
                dampg_vect=[];

                % Calculate the number of unique pair combinations
                n = length(l_vect); % Number of elements in the vector
                k = 2; % Choosing two elements
                unique_combs = nchoosek(n, k);

                % Add the repetitions (each element paired with itself)
                repetitions = n;

                % Total combinations
                bld_cases = unique_combs + repetitions;

            elseif strcmp(rf_fldr, 'Bld_with_Walls')
                l_vect=5;
                b_vect=5;
                h = 3;
                wall_config = [1,2,3,4,5,6,7,8,9,10];
                % wall_config = [1,2];
                dampg_vect=[];
                bld_cases = length(wall_config);
            elseif strcmp(rf_fldr, 'Vary_DampRatio')
                l_vect=5;
                b_vect=5;
                h = 3;
                wall_config=[];
                dampg_vect = 0:0.005:0.1;
                bld_cases=length(dampg_vect);
            end
        end
        %%
        function [B_f,L_f]=get_LfBf(ftyp,n_esize)
            if strcmp(ftyp,'PLATE')
                B_f = n_esize/2;
                L_f = n_esize/2;
            else
                B_f = 0.75;
                L_f = 0.75;
            end
        end
        %%
        function [n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=get_nstr_nrxy_fndn_soil_info(bld_soil_fndn)
            if strcmp(bld_soil_fndn, 'nstr3_Plate450')
                n_str = 3;
                n_rx = 2;
                n_ry = 3;
                ftyp = 'PLATE';
                V_s = 450;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            elseif strcmp(bld_soil_fndn, 'nstr2_Plate450')
                n_str = 2;
                n_rx = 2;
                n_ry = 3;
                ftyp = 'PLATE';
                V_s = 450;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            elseif strcmp(bld_soil_fndn, 'nstr2_Plate100')
                n_str = 2;
                n_rx = 2;
                n_ry = 3;
                ftyp = 'PLATE';
                V_s = 100;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            elseif strcmp(bld_soil_fndn, 'nstr1_Footing450')
                n_str = 1;
                n_rx = 1;
                n_ry = 1;
                ftyp = 'FOOTING';
                V_s = 450;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            elseif strcmp(bld_soil_fndn, 'nstr1_Footing100')
                n_str = 1;
                n_rx = 1;
                n_ry = 1;
                ftyp = 'FOOTING';
                V_s = 100;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            elseif strcmp(bld_soil_fndn, 'nstr3_Plate100')
                n_str = 3;
                n_rx = 2;
                n_ry = 3;
                ftyp = 'PLATE';
                V_s = 100;
                n_esize=0.5;
                [B_f,L_f]=fns_Inpt_BldPara.get_LfBf(ftyp,n_esize);
            end
        end

    end
end