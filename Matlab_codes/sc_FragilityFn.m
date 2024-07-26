%% This code performs an analysis on seismic data across 
% multiple building configurations and event scenarios.
% It retrieves event information and preallocates arrays 
% for storing results. For each building configuration, 
% it fetches relevant parameters and iterates over 
% different event scenarios. 
% For each event, it retrieves ground motion data 
% and calculates peak ground velocities (PGVs) and 
% maximum response values in multiple directions. 
% These values are aggregated and analyzed to identify 
% exceedances of predefined limits (KB values from DIN 4150-2 
% and v_max from DIN 4150-3). 
% The code then sorts the PGV data, 
% checks for exceedances against horizontal and vertical limits,
% and saves the results into a MATLAB file for further use.
%%
clear; clc; close all;  
rf_fldr_vect = {'MultiUnitBld_GeomVary_3lby2', 'Bld_with_Walls', 'Vary_DampRatio'};  % Define building configuration folders
data_set = 'Insheim_1';  % Define data set
%----------------------------------------------------%
[event_vect, nzero, bf_nm_u, bf_nm_v, cols, s_dir, n_snr, cmpt] = ...
    fns_EvntData.getEventList(data_set);  % Retrieve event list and related information
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
%%
for irf=1:length(rf_fldr_vect)  % Loop over each building configuration
    rf_fldr=rf_fldr_vect{irf};  % Get current building configuration
    disp(['selectedRfFldr: ', rf_fldr])
    [l_vect,b_vect,h,wall_config,dampg_vect,bldcases]=...
        fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);  % Get building parameters

    if strcmp(rf_fldr, 'MultiUnitBld_GeomVary_3lby2')  % Check specific configuration
        bld_soil_fndn_vect ={'nstr2_Plate450','nstr3_Plate450'};
    else
        bld_soil_fndn_vect ={'nstr3_Plate450'};
    end

    for ibld=1:length(bld_soil_fndn_vect)  % Loop over each soil-foundation type
        bld_soil_fndn=bld_soil_fndn_vect{ibld};  % Get current soil-foundation type
        disp(['bld_soil_fndnPara: ', bld_soil_fndn])
        [n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=...
            fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);  % Get soil-foundation parameters

        % Initialize arrays for event-specific data
        all_PGVx = [];
        all_PGVy = [];
        all_PGVz = [];
        all_max_Vx = [];
        all_max_Vy = [];
        all_max_Vz = [];
        all_max_KBx = [];
        all_max_KBy = [];
        all_max_KBz = [];

        for ie = 1:length(event_vect)  % Loop over each event
            evnt = event_vect{ie};  % Get current event
            disp(['evnt: ', evnt]);
            [stn_vect, date_evnt, time_evnt, r_vect] =...
                fns_data_process.get_event_fordataprocess(data_set, evnt);  % Get event data for processing
            n_stns = length(stn_vect);  % Number of stations

            % Preallocate arrays for PGV and max response values
            PGVx_vect = zeros(n_stns,1);
            PGVy_vect = zeros(n_stns,1);
            PGVz_vect = zeros(n_stns,1);
            max_Vx = zeros(n_stns,bldcases);
            max_Vy = zeros(n_stns,bldcases);
            max_Vz = zeros(n_stns,bldcases);
            max_KBx = zeros(n_stns,bldcases);
            max_KBy = zeros(n_stns,bldcases);
            max_KBz = zeros(n_stns,bldcases);

            [f_x_max, f_y_max, f_z_max, max_Vxmat, max_Vymat,...
                max_Vzmat,t_x_max, t_y_max, t_z_max,...
                max_Vx_KB_f_mat, max_Vy_KB_f_mat,max_Vz_KB_f_mat] =...
                fns_CloudAnalysis.import_DINvals(data_set,...
                date_evnt, time_evnt, rf_fldr, bld_soil_fndn);  % Import cloud analysis data

            for i_stn = 1:n_stns  % Loop over each station
                stn = stn_vect{i_stn};  % Get current station
                [ff_fldr] = fns_EvntData.get_ff_fldr1(data_set, stn,...
                    date_evnt, time_evnt);  % Get folder for station data
                [PGVx, PGVy, PGVz] = fns_CloudAnalysis.import_PGV(ff_fldr);  % Import PGV values
                PGVx_vect(i_stn) = PGVx;
                PGVy_vect(i_stn) = PGVy;
                PGVz_vect(i_stn) = PGVz;
                % Extract max values for the current station and building configuration
                max_Vx(i_stn,:) = max_Vxmat(i_stn,:,n_str+1);
                max_Vy(i_stn,:) = max_Vymat(i_stn,:,n_str+1);
                max_Vz(i_stn,:) = max_Vzmat(i_stn,:,n_str+1);
                max_KBx(i_stn,:) = max_Vx_KB_f_mat(i_stn,:,n_str+1);
                max_KBy(i_stn,:) = max_Vy_KB_f_mat(i_stn,:,n_str+1);
                max_KBz(i_stn,:) = max_Vz_KB_f_mat(i_stn,:,n_str+1);
            end

            % Aggregate PGV and max response values for all stations and events
            all_PGVx = [all_PGVx; PGVx_vect];
            all_PGVy = [all_PGVy; PGVy_vect];
            all_PGVz = [all_PGVz; PGVz_vect];
            all_max_Vx = [all_max_Vx; max_Vx];
            all_max_Vy = [all_max_Vy; max_Vy];
            all_max_Vz = [all_max_Vz; max_Vz];
            all_max_KBx = [all_max_KBx; max_KBx];
            all_max_KBy = [all_max_KBy; max_KBy];
            all_max_KBz = [all_max_KBz; max_KBz];
        end

        % Store aggregated max response values for each building configuration
        evry_max_Vx = [evry_max_Vx all_max_Vx];
        evry_max_Vy = [evry_max_Vy all_max_Vy];
        evry_max_Vz = [evry_max_Vz all_max_Vz];
        evry_max_KBx = [evry_max_KBx all_max_KBx];
        evry_max_KBy = [evry_max_KBy all_max_KBy];
        evry_max_KBz = [evry_max_KBz all_max_KBz];
    end
end

% Concatenate the matrices along the third dimension
KB3D_mat = cat(3, evry_max_KBx, evry_max_KBy, evry_max_KBz);
% Find the maximum along the third dimension
KBmax_mat = max(KB3D_mat, [], 3);
%%
nGM = size(KBmax_mat, 1);  % Number of ground motions
nbld = size(KBmax_mat, 2);  % Number of buildings
limKB_vect = [0.15 3 0.1 0.2];  % Define limit values for KB

% Initialize exceedance matrix
KB_exceed_mat = zeros(length(limKB_vect), nGM);
for i = 1:length(limKB_vect)  % Loop over each limit value
    limKB = limKB_vect(i);  % Current limit value
    fraction_lowlimKB_x = zeros(1, nGM);
    % Calculate the fraction of buildings exceeding the limit for each ground motion
    for j = 1:nGM
        fraction_lowlimKB_x(j) = sum(KBmax_mat(j, :) >= limKB);
        KB_exceed_mat(i, j) = fraction_lowlimKB_x(j);
    end
end
%%
% Sort PGV values
sorted_PGVx = (sort(all_PGVx)).';
sorted_PGVy = (sort(all_PGVx)).';
sorted_PGVz = (sort(all_PGVx)).';
sorted_PGV = [sorted_PGVx; sorted_PGVy; sorted_PGVz];
%%
lim_vhor = 15e-3;  % Horizontal velocity limit
lim_vvert = 20e-3;  % Vertical velocity limit
vx_exceed_mat = zeros(1, nGM);  % Initialize exceedance matrices
vy_exceed_mat = zeros(1, nGM);
vz_exceed_mat = zeros(1, nGM);

for j = 1:nGM  % Loop over each ground motion
    vx_exceed_mat(j) = sum(evry_max_Vx(j, :) >= lim_vhor);  % Count exceedances for horizontal velocity
    vy_exceed_mat(j) = sum(evry_max_Vy(j, :) >= lim_vhor);  % Count exceedances for horizontal velocity
    vz_exceed_mat(j) = sum(evry_max_Vz(j, :) >= lim_vvert);  % Count exceedances for vertical velocity
end

v_exceed_mat = [vx_exceed_mat; vy_exceed_mat; vz_exceed_mat];  % Combine exceedance matrices

%% Save the results
cd SAVE_DATA
cd Fragility_Fn
save fragilitydata.mat sorted_PGV KB_exceed_mat v_exceed_mat nbld  % Save data to file
cd ..
cd ..
