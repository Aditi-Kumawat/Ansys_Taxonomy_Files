classdef fns_data_process
    methods (Static)
        %%
        function [fn_fft_ss,freq_ss,fn_ifft]= fun_fftandifft(t_in,Fs,fn_in)
            df=(Fs/length(t_in));
            fn_fft = fft(fn_in) * (1/Fs);
            n_fft=length(fn_fft);
            fn_fft_ss =fn_fft(1:n_fft/2+1,:);
            freq = 0:df:Fs;
            freq_ss = freq(1:n_fft/2+1);
            fn_fft_pad = [fn_fft_ss; conj(flipud(fn_fft_ss(2:end-1,:)))];
            fn_ifft = ifft(fn_fft_pad*Fs, n_fft ,1, 'symmetric');
        end
        %%
        function [fn_fft_ss,freq,fn_ifft,t]= fun_fftandifft_old(t_in,Fs,fn_in)
            nfft = 2^nextpow2(length(t_in));
            freq = Fs / 2 * linspace(0, 1, nfft/2+1);
            fn_fft = fft(fn_in) * (1/Fs);
            %             fn_fft_ss = 2 * fn_fft(1:nfft/2+1,:); %% orignal
            fn_fft_ss =fn_fft(1:nfft/2+1,:);
            %     u_fft_ss=v_fft_ss./1i./(2*pi*f_matrix);

            fn_fft_pad = [fn_fft_ss; conj(flipud(fn_fft_ss(2:end-1,:)))];
            %             fn_ifft = (0.5)*ifft(fn_fft_pad*Fs, nfft, 1, 'symmetric'); %
            fn_ifft = ifft(fn_fft_pad*Fs, nfft, 1, 'symmetric');
            fn_ifft = fn_ifft(1:length(t_in),:); % Truncate to same length as v
            t = (0:length(fn_ifft)-1) / Fs;
        end

        %%
        function save_data_time(qnt, stn, date, time_evnt,ff_fldrnew,t,fn_time)
            for i_v=1:3
                cd ..
                file_name_2 = strcat(sprintf('%s_%d_%s_%s_%s',qnt, i_v, stn, date, time_evnt), '.txt');
                full_path = fullfile(ff_fldrnew, file_name_2);
                disp(['Trying to open: ', full_path]);
                fileID = fopen(full_path, 'w');
                fprintf(fileID, 'Transient signals (denoised)  \n');
                fprintf(fileID, 'Units as in the International System of Units (Hz, m, m/s, m/s^2) \n');
                fprintf(fileID, 'Station, UTC, component values in the filename \n');
                fprintf(fileID, '-------------------- \n');

                fprintf(fileID, 't(s)\tRe\n');

                for i_t = 1:length(t)
                    fprintf(fileID, ' %14.7e %14.7e\n', t(i_t), real(fn_time(i_t,i_v)));
                end

                fclose(fileID);
                cd Matlab_codes
            end
        end
        function save_data_freq(qnt, stn, date, time_evnt,ff_fldrnew,freq,fn_fft_ss)
            for i_v=1:3
                cd ..
                file_name_2 = strcat(sprintf('fft%s_%d_%s_%s_%s',qnt, i_v, stn, date, time_evnt), '.txt');
                full_path = fullfile(ff_fldrnew, file_name_2);
                disp(['Trying to open: ', full_path]);
                fileID = fopen(full_path, 'w');
                fprintf(fileID, 'frequency Spectra (only positive half of the frequency range, not doubled in amplitude) \n');
                fprintf(fileID, 'Units as in the International System of Units (Hz, m, m/s, m/s^2) \n');
                fprintf(fileID, 'Station, UTC, component values in the filename \n');
                fprintf(fileID, '-------------------- \n');

                fprintf(fileID, 'f(Hz)\tRe\tIm\tAbs\n');

                for i_f = 1:length(freq)
                    fprintf(fileID, ' %14.7e %14.7e %14.7e %14.7e\n', freq(i_f), real(fn_fft_ss(i_f,i_v)), imag(fn_fft_ss(i_f,i_v)), abs(fn_fft_ss(i_f,i_v)));
                end

                fclose(fileID);
                cd Matlab_codes
            end
        end

        function [stn_vect,date,time_evnt,r_vect]=get_event_fordataprocess(evnt)
            %%
            if strcmp(evnt, '2009')
                stn_vect={'LDAU'};
                date='2009_10_18';
                time_evnt='19_12_12';
                r_vect=4.66;
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2010_1')
                stn_vect={'LDAU'};
                date='2010_04_07';
                time_evnt='09_04_04';
                r_vect=5.27;
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2010_2')
                stn_vect={'LDAU'};
                date='2010_04_07';
                time_evnt='13_46_21';
                r_vect=5.27;
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2012_1')
                stn_vect={'INS2','INS3','INS4','INS5','INS6', 'INSH'};
                date='2012_11_12';
                time_evnt='11_15_04';
                r_vect=[3.19 7.67 6.76 3.52 6.67 0.64];
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2012_2')
                stn_vect={'INS2','INS3','INS4','INS5','INS6', 'INSH'};
                date='2012_11_12';
                time_evnt='12_53_02';
                r_vect=[2.88 7.36 6.84 3.63 6.57 0.95];
                % tmin=0;
                % tmax=50;

            elseif strcmp(evnt, '2013Jan')
                % stn_vect={'INS1'};
                stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7', 'INSH'};
                date='2013_01_26';
                time_evnt='19_48_27';
                r_vect=[7.86 4.55 8.65 5.00 4.96 8.54 4.3 1.6];
                % tmin=0;
                % tmax=60;

            elseif strcmp(evnt, '2013Feb') %%
                stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7','INSH'};
                date='2013_02_17';
                time_evnt='20_07_15';
                r_vect=[7.82 4.74 8.79 4.78 5.16 8.77 4.36 1.82];
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2013Oct') %% higest magnitude
                % stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
                %     'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
                stn_vect={'INSH','TMO54','INS5'};
                r_vect=[1.9 3 5.2];
                date='2013_10_02';
                time_evnt='01_13_26';
                % r_vect=[7.84 4.85 8.88 4.7 5.24 8.88 4.37 1.9 7.47 6.8 3 7.9];
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2013Nov18')
                stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
                    'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
                date='2013_11_18';
                time_evnt='12_54_15';
                r_vect=[7.27 2.8 7.0 6.2 4.5 7.2 5.1 1.49 7.98 8.82 1.35 9.66];
                % tmin=0;
                % tmax=60;

            elseif strcmp(evnt, '2013Nov21') %%
                stn_vect={'INS1','INS2','INS3','INS4','INS5','INS6','INS7',...
                    'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55'};
                date='2013_11_21';
                time_evnt='14_15_22';
                r_vect=[7.84 3.64 7.92 5.87 4.3 7.58 4.36 1.49 8.14 8.0 2.17 9.20];
                % tmin=0;
                % tmax=40;

            elseif strcmp(evnt, '2016_1')
                stn_vect={'A127A','INS3','INS4B','INS5','INS6B','INS7','INS8',...
                    'INSH', 'TMO20', 'TMO22', 'TMO54', 'TMO55', 'TMO57', 'TMO66'};
                date='2016_02_12';
                time_evnt='06_26_04';
                r_vect=[7.40 5.91 4.32 7.48 4.77 13.11 1.16 8.11 8.52 1.71 9.53 9.73 11.0 13.71];
                % tmin=0;
                % tmax=50;

            elseif strcmp(evnt, '2016_2')
                stn_vect={'INS1','INS3','INS4B','INS5','INS6B','INS7','INS8',...
                    'INSH', 'TMO20', 'TMO54', 'TMO55', 'TMO57', 'TMO66'};
                date='2016_07_14';
                time_evnt='17_49_10';
                r_vect=[7.99 8.02 5.76 4.17 7.71 4.22 12.74 0.77 8.27 2.30 9.25 9.18 10.98];
                % tmin=0;
                % tmax=40;
            end
        end
    end
end