clear; clc; close all;
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
        all_max_Vx = [];
        all_max_Vy = [];
        all_max_Vz = [];
        all_max_KBx = [];
        all_max_KBy = [];
        all_max_KBz = [];
        for ie = 1:length(event_vect)
            evnt = event_vect{ie};
            disp(['evnt: ', evnt]);
            [stn_vect, date_evnt, time_evnt, r_vect] =...
                fns_data_process.get_event_fordataprocess(evnt);
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

ha_cl=@colors;
lcol = {ha_cl('cinnamon'),ha_cl('denim'),...
                ha_cl('crimson'),ha_cl('gray'),ha_cl('black')};
% col=ha_col('debian red','denim','dark sea green','amethyst','black');
limKB_vect=[3 0.15 0.2 0.1];
lim_vhor=15e-3;
lim_vvert=20e-3;
ndeg=3;
legend_vect={'$KB_{U}$ (day)','$KB_{L}$ (day)', '$KB_{U}$ (night)',...
    '$KB_{L}$ (night)', '$v_{max}$'};
% legend_vect={'$KB_{U}$ (day)', '$KB_{U}$ (day)','$KB_{L}$ (day)',...
%     '$KB_{L}$ (day)', '$KB_{U}$ (night)', '$KB_{U}$ (night)',...
%     '$KB_{L}$ (night)', '$KB_{L}$ (night)', '$v_{max}$', '$v_{max}$'};
fns_FragilityFns.plot_direction('x',limKB_vect,evry_max_KBx,...
    all_PGVx,ndeg,evry_max_Vx,lim_vhor,lcol,legend_vect)
fns_FragilityFns.plot_direction('y',limKB_vect,evry_max_KBy,...
    all_PGVy,ndeg,evry_max_Vy,lim_vhor,lcol,legend_vect)
fns_FragilityFns.plot_direction('z',limKB_vect,evry_max_KBz,...
    all_PGVz,ndeg,evry_max_Vz,lim_vvert,lcol,legend_vect)

