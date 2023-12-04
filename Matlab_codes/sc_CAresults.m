clear; clc; close all;

%----------------------------------------------------%
rf_fldr = fns_Inpt_BldPara.selectRfFldr();
disp(['selectedRfFldr: ', rf_fldr])
%----------------------------------------------------%
bld_soil_fndn = fns_Inpt_BldPara.select_bldsoilpara();
disp(['bld_soil_fndnPara: ', bld_soil_fndn])
%----------------------------------------------------%
data_set = fns_EvntData.select_event_stn();
disp(['selected_dataset: ', data_set])
%----------------------------------------------------%
[l_vect,b_vect,h,wall_config,dampg_vect,bldcases]=...
    fns_Inpt_BldPara.get_lbh_bldcases_for_rf_fldr(rf_fldr);
[n_str,n_rx,n_ry,V_s,ftyp,B_f,L_f]=...
    fns_Inpt_BldPara.get_nstr_nrxy_fndn_soil_info(bld_soil_fndn);
[event_vect, nzero, bf_nm_u, bf_nm_v, cols, s_dir, n_snr, cmpt] = ...
    fns_EvntData.getEventList(data_set);
%----------------------------------------------------%
% Preallocate arrays for plotting
all_PGVx = [];
all_max_Vx = [];
all_PGVy = [];
all_max_Vy = [];
all_PGVz = [];
all_max_Vz = [];
event_labels = [];
bld_labels= [];

for ie = 1:length(event_vect)
    evnt = event_vect{ie};
    disp(['evnt: ', evnt]);
    [stn_vect, date_evnt, time_evnt, r_vect] = fns_data_process.get_event_fordataprocess(evnt);
    n_stns = length(stn_vect);
    PGVx_vect = zeros(n_stns,1);PGVy_vect = zeros(n_stns,1);PGVz_vect = zeros(n_stns,1);
    [f_x_max, f_y_max, f_z_max, max_Vxmat, max_Vymat, max_Vzmat, ...
        t_x_max, t_y_max, t_z_max, max_Vx_KB_f_mat, max_Vy_KB_f_mat,...
        max_Vz_KB_f_mat] = fns_CloudAnalysis.import_DINvals(data_set,...
        date_evnt, time_evnt, rf_fldr, bld_soil_fndn);

    for i_stn = 1:n_stns
        stn = stn_vect{i_stn};
        [ff_fldr] = fns_EvntData.get_ff_fldr1(data_set, stn, date_evnt, time_evnt);
        [PGVx, PGVy, PGVz] = fns_CloudAnalysis.import_PGV(ff_fldr);
        PGVx_vect(i_stn) = PGVx;
        PGVy_vect(i_stn) = PGVy;
        PGVz_vect(i_stn) = PGVz;
        for ibld = 1:bldcases
            for istr = n_str+1

                max_Vx = max_Vxmat(i_stn,ibld,istr);
                max_Vy = max_Vymat(i_stn,ibld,istr);
                max_Vz = max_Vzmat(i_stn,ibld,istr);
            end
            % Store values for plotting
            all_PGVx = [all_PGVx; PGVx_vect(i_stn)];
            all_max_Vx = [all_max_Vx; max_Vx];
            all_PGVy = [all_PGVy; PGVy_vect(i_stn)];
            all_max_Vy = [all_max_Vy; max_Vy];
            all_PGVz = [all_PGVz; PGVz_vect(i_stn)];
            all_max_Vz = [all_max_Vz; max_Vz];
            event_labels = [event_labels; ie]; % Store event index for coloring
            bld_labels = [bld_labels; ibld]; % Store event index for coloring
        end
    end
end
% Create a custom legend
% Plotting
evntlbl_length = jet(length(event_vect)); 
bldlbl_length = jet((bldcases)); 
% Plot for PGVx vs Max Vx
figure;
scatter(all_PGVx, all_max_Vx, [], bld_labels); % Your original scatter plot
% colormap(bldlbl_length); % Define a colormap according to event length
%----------------------------------------------------%
figure;
scatter(all_PGVx, all_max_Vx, [], event_labels);
colormap(evntlbl_length);
% Create a custom legend
hold on;
hold on;
legend_entries = gobjects(length(event_vect), 1); % Initialize an array of graphic objects
for ie = 1:length(event_vect)
    legend_entries(ie) = scatter(NaN, NaN, [], evntlbl_length(ie, :), 'DisplayName', event_vect{ie});
end
legend(legend_entries, event_vect);
hold off;
%----------------------------------------------------%
figure;
scatter(all_PGVy, all_max_Vy, [], event_labels);
colormap(evntlbl_length);
% Create a custom legend
hold on;
legend_entries = gobjects(length(event_vect), 1); % Initialize an array of graphic objects
for ie = 1:length(event_vect)
    legend_entries(ie) = scatter(NaN, NaN, [], evntlbl_length(ie, :), 'DisplayName', event_vect{ie});
end
legend(legend_entries, event_vect);
hold off;
%----------------------------------------------------%
figure;
scatter(all_PGVz, all_max_Vz, [], event_labels);
colormap(evntlbl_length);
% Create a custom legend
hold on;
legend_entries = gobjects(length(event_vect), 1); % Initialize an array of graphic objects
for ie = 1:length(event_vect)
    legend_entries(ie) = scatter(NaN, NaN, [], evntlbl_length(ie, :), 'DisplayName', event_vect{ie});
end
legend(legend_entries, event_vect);
hold off;




