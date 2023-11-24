%% fns_Inpt_Para
classdef fns_EvntData
    methods (Static)
        function selectedOption = select_event_stn()
            options = {'Insheim', 'Poing', 'Unterhaching'};
            [choice, isOk] = listdlg('PromptString', 'Select an event:', ...
                'SelectionMode', 'single', ...
                'ListString', options,'ListSize',[150,150]);

            if isOk
                selectedOption = options{choice};
            else
                selectedOption = ''; % or handle the case where no selection is made
            end
        end

        %%
        function [evnt,stn,stn_vect,date,time,nzero,ff_fldr,...
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
                date='2016_12_20';
                time='03_30_51';
                nzero=6;
                ff_fldr = fullfile('GM','GM_POI2016',stn);
            elseif strcmp(name_evnt, 'Unterhaching')
                evnt='Part1';
                stn='UH1';
                stn_vect={'UH1', 'UH2', 'UH3'};
                date='2013_04_16';
                time='21_51_42';
                nzero=4;
                fldr_nm = [stn, '_', evnt];
                ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
            elseif strcmp(name_evnt, 'Insheim')
                evnt='Insheim_Oct2013';
                stn='TMO54';
                stn_vect={'INSH','TMO54','INS5'};
                date='2013_10_02';
                time='01_13_26';
                nzero=6;
                ff_fldr = fullfile('GM','GM_Insheim_Oct2013', stn);
            end
        end
        %%
        function [ff_fldr]=get_ff_fldr(name_evnt,stn)
            if strcmp(name_evnt, 'Poing')
                ff_fldr = fullfile('GM','GM_POI2016',stn);
            elseif strcmp(name_evnt, 'Unterhaching')
                fldr_nm = [stn, '_', evnt];
                ff_fldr = fullfile('GM', 'GM_UH',fldr_nm);
            elseif strcmp(name_evnt, 'Insheim')
                ff_fldr = fullfile('GM','GM_Insheim_Oct2013', stn);
            end
        end
        
    end
end