clear; clc; close all;
%% Plots the cloud analysis results for all cases of building models
% This plots added to the paper has analyzed all possible cases for the
% building model amounting to 87 building cases. The plots show the
% maximum values of the velocity and KB levels at selected locations of the
% models with the PGVs of the GMs of the selected dataset of input seismic
% data. 
% !!The results have been saved by sc_CA.m!!

rf_fldr_vect = {'MultiUnitBld_GeomVary_3lby2', 'Bld_with_Walls',...
    'Vary_DampRatio'};
data_set = 'Insheim_1';
%----------------------------------------------------%
[event_vect, nzero, bf_nm_u, bf_nm_v, cols, s_dir, n_snr, cmpt] = ...
    fns_EvntData.getEventList(data_set);
%----------------------------------------------------%
% Preallocate arrays for plotting
event_labels = [];
bld_labels= [];
evry_max_Vx = [];
evry_max_Vy = [];
evry_max_Vz = [];
evry_max_KBx = [];
evry_max_KBy = [];
evry_max_KBz = [];
for irf=1:length(rf_fldr_vect)
    rf_fldr=rf_fldr_vect{irf};
    disp(['selectedRfFldr: ', rf_fldr])
    [l_vect,b_vect,h,wall_config,dampg_vect,bldcases]=...
        fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);
    if strcmp(rf_fldr, 'MultiUnitBld_GeomVary_3lby2')
        bld_soil_fndn_vect ={'nstr2_Plate450','nstr3_Plate450'};
    else
        bld_soil_fndn_vect ={'nstr3_Plate450'};
    end
    for ibld=1:length(bld_soil_fndn_vect)
        bld_soil_fndn=bld_soil_fndn_vect{ibld};
        disp(['bld_soil_fndnPara: ', bld_soil_fndn])
        [n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=...
            fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);
        all_PGVx = [];
        all_PGVy = [];
        all_PGVz = [];
        event_labels_full = [];
        all_max_Vx = [];
        all_max_Vy = [];
        all_max_Vz = [];
        all_max_KBx = [];
        all_max_KBy = [];
        all_max_KBz = [];
        event_labels = [];
        for ie = 1:length(event_vect)
            evnt = event_vect{ie};
            disp(['evnt: ', evnt]);
            [stn_vect, date_evnt, time_evnt, r_vect] =...
                fns_data_process.get_event_fordataprocess(data_set,evnt);
            n_stns = length(stn_vect);
            PGVx_vect = zeros(n_stns,1);
            PGVy_vect = zeros(n_stns,1);
            PGVz_vect = zeros(n_stns,1);
            max_Vx = zeros(n_stns,bldcases);
            max_Vy = zeros(n_stns,bldcases);
            max_Vz = zeros(n_stns,bldcases);
            max_KBx = zeros(n_stns,bldcases);
            max_KBy = zeros(n_stns,bldcases);
            max_KBz = zeros(n_stns,bldcases);
            ie_vect=zeros(n_stns,1);
            [f_x_max, f_y_max, f_z_max, max_Vxmat, max_Vymat,...
                max_Vzmat,t_x_max, t_y_max, t_z_max,...
                max_Vx_KB_f_mat, max_Vy_KB_f_mat,max_Vz_KB_f_mat] =...
                fns_CloudAnalysis.import_DINvals(data_set,...
                date_evnt, time_evnt, rf_fldr, bld_soil_fndn);

            for i_stn = 1:n_stns
                stn = stn_vect{i_stn};
                [ff_fldr] = fns_EvntData.get_ff_fldr1(data_set, stn,...
                    date_evnt, time_evnt);
                [PGVx, PGVy, PGVz] = fns_CloudAnalysis.import_PGV(ff_fldr);
                PGVx_vect(i_stn) = PGVx;
                PGVy_vect(i_stn) = PGVy;
                PGVz_vect(i_stn) = PGVz;
                ie_vect(i_stn) = ie;
                max_Vx(i_stn,:) = max_Vxmat(i_stn,:,n_str+1);
                max_Vy(i_stn,:) = max_Vymat(i_stn,:,n_str+1);
                max_Vz(i_stn,:) = max_Vzmat(i_stn,:,n_str+1);
                max_KBx(i_stn,:) = max_Vx_KB_f_mat(i_stn,:,n_str+1);
                max_KBy(i_stn,:) = max_Vy_KB_f_mat(i_stn,:,n_str+1);
                max_KBz(i_stn,:) = max_Vz_KB_f_mat(i_stn,:,n_str+1);
            end
            all_PGVx = [all_PGVx; PGVx_vect];
            all_PGVy = [all_PGVy; PGVy_vect];
            all_PGVz = [all_PGVz; PGVz_vect];
            event_labels_full=[event_labels_full; ie_vect];
            all_max_Vx = [all_max_Vx; max_Vx];
            all_max_Vy = [all_max_Vy; max_Vy];
            all_max_Vz = [all_max_Vz; max_Vz];
            all_max_KBx = [all_max_KBx; max_KBx];
            all_max_KBy = [all_max_KBy; max_KBy];
            all_max_KBz = [all_max_KBz; max_KBz];
        end
        evry_max_Vx = [evry_max_Vx all_max_Vx];
        evry_max_Vy = [evry_max_Vy all_max_Vy];
        evry_max_Vz = [evry_max_Vz all_max_Vz];
        evry_max_KBx = [evry_max_KBx all_max_KBx];
        evry_max_KBy = [evry_max_KBy all_max_KBy];
        evry_max_KBz = [evry_max_KBz all_max_KBz];
    end
end

ylbl_vect1={'$\max[v_x(t)]$,~mm/s', '$\max[v_y(t)]$,~mm/s',...
    '$\max[v_z(t)]$,~mm/s'};
% ylbl_vect2={'$KB_{Fmax}(x)$', '$KB_{Fmax}(y)$', '$KB_{Fmax}(z)$'};
ylbl_vect2={'$\max[KB_{F}]$','$\max[KB_{F}]$','$\max[KB_{F}]$'};
% xlbl_vect={'PGV(x),~m/s', 'PGV(y),~m/s', 'PGV(z),~m/s'};
xlbl_vect = {'log(PGV($x$)), mm/s', 'log(PGV($y$)), mm/s',...
    'log(PGV($z$)), mm/s'};
filename = {'VmaxPlot_x.pdf','VmaxPlot_y.pdf','VmaxPlot_z.pdf',...
    'KBmaxPlot_x.pdf','KBmaxPlot_y.pdf','KBmaxPlot_z.pdf'};
filename = strcat(data_set, '_', filename);

if strcmp(data_set, 'Poing')
    event_vect_lbl = { '$2016$', '$2017$' };
    ylimits_set = [0 8; 0 3.5; 0 26; 0 4; 0 2; 0 6];
elseif strcmp(data_set, 'Insheim_1')
    event_vect_lbl = {'$2009$', '$2010_1$', '$2010_2$', '$2012_1$',...
        '$2012_2$', '$2013_1$', '$2013_2$', '$2013_3$', '$2013_4$',...
        '$2013_5$', '$2016_1$', '$2016_2$'};
    ylimits_set = [0 12; 0 6; 0 25; 0 14; 0 14; 0 14];
end
% Concatenate the matrices along the third dimension
KB3D_mat = cat(3, evry_max_KBx, evry_max_KBy, evry_max_KBz);

% Find the maximum along the third dimension
KBmax_mat = max(KB3D_mat, [], 3);
% data_vect = {evry_max_Vx, evry_max_Vy, evry_max_Vz, evry_max_KBx, evry_max_KBy, evry_max_KBz};
data_vect = {evry_max_Vx, evry_max_Vy, evry_max_Vz, KBmax_mat, KBmax_mat, KBmax_mat};
PGV_vect={all_PGVx, all_PGVy, all_PGVz, all_PGVx, all_PGVy, all_PGVz};
fct_mm_set = [1e3, 1e3, 1e3, 1, 1, 1];

for i = 1:6
    ylimits = ylimits_set(i, :);
    fct_mm = fct_mm_set(i);

    if i <= 3
        ylbl_vect = ylbl_vect1;
    else
        ylbl_vect = ylbl_vect2;
    end

    fns_CloudAnalysis.plot_CAfull(event_vect_lbl, PGV_vect{i}*1e3, data_vect{i} * fct_mm, ...
        event_labels_full, ylbl_vect{mod(i-1,3)+1}, xlbl_vect{mod(i-1,3)+1}, ...
        filename{i}, ylimits);
end


