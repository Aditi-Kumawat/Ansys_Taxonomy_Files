%% fns_Inpt_Para
classdef fns_EvntData
    methods (Static)
        function selectedOption = select_event_stn()
            options = {'Insheim_1', 'Poing', 'Unterhaching'};
            [choice, isOk] = listdlg('PromptString', 'Select a location:', ...
                'SelectionMode', 'single', ...
                'ListString', options,'ListSize',[150,150]);

            if isOk
                selectedOption = options{choice};
            else
                selectedOption = ''; % or handle the case where no selection is made
            end
        end

        %%
        function [evnt,stn,stn_vect,date_evnt,time,nzero,ff_fldr,...
                bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=get_event_stn(name_evnt)
            bf_nm_u = 'fftd_%d_%s_%s_%s';
            bf_nm_v = 'fftv_%d_%s_%s_%s';
            cols = {'Freq', 'Re','Im','Amp'};
            s_dir = [1 2 3];
            n_snr = numel(s_dir);
            cmpt = {'X', 'Y', 'Z'};
            if strcmp(name_evnt, 'Poing')
                evnt='Po2016';
                stn='POI01';
                stn_vect={'POI01', 'POI02', 'POI03'};
                date_evnt='2016_12_20';
                time='03_30_51';
                nzero=6;
                ff_fldr = fullfile('GM','GM_POI2016',stn);
            elseif strcmp(name_evnt, 'Unterhaching')
                evnt='Part1';
                stn='UH1';
                stn_vect={'UH1', 'UH2', 'UH3'};
                date_evnt='2013_04_16';
                time='21_51_42';
                nzero=4;
                fldr_nm = [stn, '_', evnt];
                ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
            elseif strcmp(name_evnt, 'Insheim_1')
                evnt='Insheim_Oct2013';
                stn='TMO54';
                stn_vect={'INSH','TMO54','INS5'};
                date_evnt='2013_10_02';
                time='01_13_26';
                nzero=6;
                ff_fldr = fullfile('GM','GM_Insheim_Oct2013', stn);
            end
        end
        %%
        function [event_vect,nzero,bf_nm_u,bf_nm_v,cols,s_dir,n_snr,cmpt]=getEventList(name_dataset)
            bf_nm_u = 'fftd_%d_%s_%s_%s';
            bf_nm_v = 'fftv_%d_%s_%s_%s';
            cols = {'Freq', 'Re','Im','Amp'};
            s_dir = [1 2 3];
            n_snr = numel(s_dir);
            cmpt = {'X', 'Y', 'Z'};
            if strcmp(name_dataset, 'Poing')

            elseif strcmp(name_dataset, 'Unterhaching')
            elseif strcmp(name_dataset, 'Insheim_1')
                event_vect = {'2009', '2010_1', '2010_2', '2012_1', '2012_2', '2013_1', '2013_2', '2013_3', '2013_4', '2013_5', '2016_1', '2016_2'};
                nzero=6;
            end
        end

        %%
        function [ff_fldr]=get_ff_fldr(name_evnt,stn)
            if strcmp(name_evnt, 'Poing')
                ff_fldr = fullfile('GM','GM_POI2016',stn);
            elseif strcmp(name_evnt, 'Unterhaching')
                fldr_nm = [stn, '_', evnt];
                ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
            elseif strcmp(name_evnt, 'Insheim_1')
                ff_fldr = fullfile('GM','GM_Insheim_Oct2013', stn);
            end
        end
        %%
        function [ff_fldr]=get_ff_fldr1(data_set,stn,date_evnt,time_evnt)
            if strcmp(data_set, 'Poing')
                ff_fldr = fullfile('GM','GM_POI2016',stn);
            elseif strcmp(data_set, 'Unterhaching')
                fldr_nm = [stn, '_', evnt];
                ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
            elseif strcmp(data_set, 'Insheim_1')
                fldr_evnt=['GM_', date_evnt, '_', time_evnt];
                fldr_stn_nm = ['GM_', stn, '_', date_evnt, '_', time_evnt];
                ff_fldr=fullfile('GM', 'GM_Insheim_1', fldr_evnt, fldr_stn_nm);
            end
        end



    end
end