clear; clc; close all;
%% Plots the cloud analysis results for individual set of results for different 
% cases of models
% The plots added to the paper: sc_CAresultsFull.m
% The results have been saved by sc_CA.m
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
all_PGVy = [];
all_PGVz = [];
all_max_Vx = [];
all_max_Vy = [];
all_max_Vz = [];
all_max_KBx = [];
all_max_KBy = [];
all_max_KBz = [];
event_labels = [];
bld_labels= [];

for ie = 1:length(event_vect)
    evnt = event_vect{ie};
    disp(['evnt: ', evnt]);
    [stn_vect, date_evnt, time_evnt, r_vect] =...
        fns_data_process.get_event_fordataprocess(data_set,evnt);
    n_stns = length(stn_vect);
    PGVx_vect = zeros(n_stns,1);
    PGVy_vect = zeros(n_stns,1);
    PGVz_vect = zeros(n_stns,1);
    [f_x_max, f_y_max, f_z_max, max_Vxmat, max_Vymat, max_Vzmat, ...
        t_x_max, t_y_max, t_z_max, max_Vx_KB_f_mat, max_Vy_KB_f_mat,...
        max_Vz_KB_f_mat] = fns_CloudAnalysis.import_DINvals(data_set,...
        date_evnt, time_evnt, rf_fldr, bld_soil_fndn);

    for i_stn = 1:n_stns
        stn = stn_vect{i_stn};
        [ff_fldr] = fns_EvntData.get_ff_fldr1(data_set, stn,...
            date_evnt, time_evnt);
        [PGVx, PGVy, PGVz] = fns_CloudAnalysis.import_PGV(ff_fldr);
        PGVx_vect(i_stn) = PGVx;
        PGVy_vect(i_stn) = PGVy;
        PGVz_vect(i_stn) = PGVz;
        for ibld = 1:bldcases
            for istr = n_str+1

                max_Vx = max_Vxmat(i_stn,ibld,istr);
                max_Vy = max_Vymat(i_stn,ibld,istr);
                max_Vz = max_Vzmat(i_stn,ibld,istr);
                max_KBx = max_Vx_KB_f_mat(i_stn,ibld,istr);
                max_KBy = max_Vy_KB_f_mat(i_stn,ibld,istr);
                max_KBz = max_Vz_KB_f_mat(i_stn,ibld,istr);
            end
            % Store values for plotting
            all_PGVx = [all_PGVx; PGVx_vect(i_stn)];
            all_PGVy = [all_PGVy; PGVy_vect(i_stn)];
            all_PGVz = [all_PGVz; PGVz_vect(i_stn)];
            all_max_Vx = [all_max_Vx; max_Vx];
            all_max_Vy = [all_max_Vy; max_Vy];
            all_max_Vz = [all_max_Vz; max_Vz];
            all_max_KBx = [all_max_KBx; max_KBx];
            all_max_KBy = [all_max_KBy; max_KBy];
            all_max_KBz = [all_max_KBz; max_KBz];
            event_labels = [event_labels; ie];
            bld_labels = [bld_labels; ibld];
        end
    end
end
bldlbl_length = jet((bldcases));
% figure;
% scatter(all_PGVx, all_max_Vx, [], bld_labels); % Your original scatter plot
% colormap(bldlbl_length); % Define a colormap according to event length
%----------------------------------------------------%
event_vect_lbl = {'$2009$', '$2010_1$', '$2010_2$', '$2012_1$',...
    '$2012_2$', '$2013_1$', '$2013_2$', '$2013_3$', '$2013_4$',...
    '$2013_5$', '$2016_1$', '$2016_2$'};

ylbl_vect1={'$v_{max}(x)$,~m/s', '$v_{max}(y)$,~m/s',...
    '$v_{max}(z)$,~m/s'};
ylbl_vect2={'$KB_{Fmax}(x)$', '$KB_{Fmax}(y)$', '$KB_{Fmax}(z)$'};
xlbl_vect={'PGV(x),~m/s', 'PGV(y),~m/s', 'PGV(z),~m/s'};
xlbl_vect = {'log(PGV($x$)), m/s', 'log(PGV($y$)), m/s',...
    'log(PGV($z$)), m/s'};
titleText = {'$v_{max}$ for Building Responses ($x$-dir)',...
    '$v_{max}$ for Building Responses ($y$-dir)',...
    '$v_{max}$ for Building Responses ($z$-dir)',...
    'KB value for Building Responses ($x$-dir)',...
    'KB value for Building Responses ($y$-dir)',...
    'KB value for Building Responses ($z$-dir)'};
filename = {'VmaxPlot_x.pdf','VmaxPlot_y.pdf','VmaxPlot_z.pdf',...
    'KBmaxPlot_x.pdf','KBmaxPlot_y.pdf','KBmaxPlot_z.pdf'};
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVx,all_max_Vx,...
    event_labels,ylbl_vect1{1},xlbl_vect{1},titleText{1},filename{1})
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVy,all_max_Vy,...
    event_labels,ylbl_vect1{2},xlbl_vect{2},titleText{2},filename{2})
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVz,all_max_Vz,...
    event_labels,ylbl_vect1{3},xlbl_vect{3},titleText{3},filename{3})
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVx,all_max_KBx,...
    event_labels,ylbl_vect2{1},xlbl_vect{1},titleText{4},filename{4})
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVy,all_max_KBy,...
    event_labels,ylbl_vect2{2},xlbl_vect{2},titleText{5},filename{5})
fns_CloudAnalysis.plot_CA(event_vect_lbl,all_PGVz,all_max_KBz,...
    event_labels,ylbl_vect2{3},xlbl_vect{3},titleText{6},filename{6})







