classdef fns_KBvalue
    methods (Static)
        %% compute KB value
        %input: 1. freq of your data
        %       2. data in frequency domain
        %       3. threshold of the filter, e.g. 5.6Hz

        function KB = fn_kb_highpass(freq, signal_freq, highpass_threshold)
            KB = signal_freq./(sqrt(1+power((highpass_threshold./freq),2)));
        end
        %%
        function KB_f = fn_rms_kb(time, KB , time_constant)
            KB_f = zeros(size(KB));
            if length(time) == length(KB)
                %initialized
                df_ = time(3) - time(2);
                %compute KB_f
                for i = 1:length(KB_f)
                    exp_term = transpose(exp(-(time(i)-time(1:i))/time_constant));
                    KB_f(i,:) = sqrt(sum(exp_term.*((KB(1:i,:).^2)*df_),1)/time_constant);
                end
            else
                disp("ERROR! different length")
            end
        end

        %% Code refer to DIN4150-2
        %% Guildline values for evaluating human exposure to vibration
        %test = find_A_value("ReinesWohngebiet","night")

        function A_values = find_A_values(location,time)
            % A_value = [Au,Ao,Ar], please check DIN 4150-2 page7 for detail

            time_ = ["day","night"];
            loc_ = ["Indistriegebiet","Gewerbegebiet","Kerngebiet","ReinesWohngebiet","Special"];
            day = {[0.4,6.0,0.2]; [0.3,6.0,0.15];[0.2,0.4,0.1]; [0.15,3.0,0.07]; [0.1,3.0,0.05]};
            night = {[0.3,0.6,0.15]; [0.2,0.4,0.1];[0.15,0.3,0.07]; [0.1,0.2,0.05]; [0.1,0.15,0.05]};
            T = table(day,night, 'RowNames', loc_, 'VariableNames', ["Day","Night"]);

            %initialize
            A_values = [0,0,0] ;

            if ismember(location,loc_)
                if ismember(time,time_)
                    if strcmp(time,time_(1)) %case = day
                        A_values = T.Day(location);
                        A_values = A_values{1};
                    else %case = night
                        A_values = T.Night(location);
                        A_values = A_values{1};
                    end
                else
                    disp("ERROR! wrong input of time varaible");
                    disp("please input:");
                    disp("  day");
                    disp("  night");
                end
            else
                disp("ERROR! wrong input of location varaible");
                disp("please input:");
                disp("    Indistriegebiet");
                disp("    Gewerbegebiet");
                disp("    Kerngebiet");
                disp("    ReinesWohngebiet");
                disp("    Special");
            end
        end

        %%
        function required_str = fn_evaluate_A_criteria(KB_f_max,KB_Ftr_,criteria,frequent_or_not)
            Au = criteria(1);
            Ao = criteria(2);
            Ar = criteria(3);
            if KB_f_max <= Au
                required_str = "Au passed";
            else
                if KB_f_max <= Ao
                    if frequent_or_not
                        if KB_Ftr_ <= Ar
                            required_str = "Ar passed";
                        else
                            required_str = "failed";
                        end
                    else
                        required_str = "Ao passed";
                    end
                else
                    required_str = "failed";
                end
            end
        end
        %%
        function KB_fmax_appr = fn_appr_evaluation_unweighted_signal(V_max,freq,cf, threshold)
            % KB_fmax_appr = [computed value, lower bound (-15%), upper bound(+15%)]
            KB_fmax_appr = zeros(length(V_max),3);
            %value = zeros(length(V_max),1);
            value = (1/sqrt(2)).*(1./sqrt(1+(power((5.6 ./freq),2))).*V_max*cf);
            KB_fmax_appr(:,1) = value;
            KB_fmax_appr(:,2) = value*(1-threshold);
            KB_fmax_appr(:,3) = value*(1+threshold);
        end
    end
end