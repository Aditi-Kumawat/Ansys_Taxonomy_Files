classdef fns_FragilityFns
    methods (Static)
        %% -------------------------------------------------------------------------------%
        function [theta, beta] = fn_mle_pc(IM, num_gms, num_collapse)
            %% the following function calculates the median (theta)
            % and lognormal standard deviation (beta) of a fragility function.
            % The method uses maximum likelihood estimation (MLE) to fit a
            % lognormal distribution to observed collapse data at various
            % intensity measure (IM) levels.
            % INPUTS:
            % IM            1xn           IM levels of interest
            % num_gms       1x1 or 1xn    number of ground motions used at each IM level
            % num_collapse 	1xn           number of collapses observed at each IM level
            %
            % OUTPUTS:
            % theta         1x1           median of fragility function
            % beta          1x1           lognormal standard deviation of fragility function

            %% Initial guess for the fragility function parameters theta and beta
            % ** Use method of moments **
            x0 = [mean(log(IM)), std(log(IM))];

            %% Run optimization
            options = optimset('MaxFunEvals',1000, 'GradObj', 'off'); %maximum 1000 iterations, gradient of the function not provided
            x = fminsearch(@mlefit, x0, options, num_gms, num_collapse, IM) ;
            theta = exp(x(1)); % return theta in linear space
            beta = x(2);

            %% Objective function to be optimized
            function [loglik] = mlefit(params, num_gms, num_collapse, IM)

                % ** Penalize any negative beta with a very large loglik value **
                if params(2) < 0
                    loglik = 1e10;

                else
                    % estimated probabilities of collapse, given the current fragility function
                    % parameter estimates
                    p = normcdf(log(IM), (params(1)), params(2));

                    % likelihood of observing num_collapse(i) collapses, given num_gms
                    % observations, using the current fragility function parameter estimates
                    likelihood = binopdf(num_collapse', num_gms', p'); %

                    % ** Cannot have zero likelihood value, so replace every zero likelihood
                    % value with the smallest positive normalized fixed-point value **
                    likelihood(likelihood == 0) = realmin;

                    % sum negative log likelihood (we take the negative value because we want
                    % the maximum log likelihood, and the function is searching for a minimum)
                    loglik = -sum(log(likelihood));
                end
            end
        end

        %% -------------------------------------------------------------------------------%
        function save_figure_to_pdf(filename)
            % save figure to PDF
            cd SAVE_FIGS
            cd Fragility_Fn
            saveas(gcf, filename);
            cd ..
            cd ..
        end
    end
end