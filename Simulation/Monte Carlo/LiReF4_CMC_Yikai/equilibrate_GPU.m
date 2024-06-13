function [init, coef, E_si, E_int] = equilibrate_GPU(const, ion, params, E_si, E_int, E_tot, hamI, coef, basis, spin_e, spin_n)
% Convert parameters and initial arrays to GPU arrays
t_lim = params.tEQ; % equilibration time limit
accpRate = zeros(params.tEQ, length(params.temp), size(params.field, 2)); % acceptance rate tracker
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
% Main equilibration loop
for tt = 1:length(params.temp)
    beta = gpuArray(1 / (const.J2meV * const.kB * params.temp(tt)));
    for ff = 1:size(params.field,2)
        field = gpuArray(params.field(:,ff));
        eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
        nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
        fprintf('Temperature: %.2f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))

        trap_check = 0; % critical slowing-down token
        % convgTkr = 0; % cluster update convergence tracker
        Et_temp = zeros(t_lim,1);
        dE = zeros(t_lim,1);
        for time = 1:t_lim
            % thermalize the electronic spins
            [E_si(:,tt,ff), E_int(:,tt,ff), dE(time), coef(:,:,tt,ff), eSpin, accpRate(time,tt,ff)] = thermalize_GPU(const, ion,...
                params, E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin, beta);
            % thermalize the nuclear spins
            % if params.hyp && mod(time,20) == 0 % wot slow down nuclear spin update
            if params.hyp
                nSpin = therm_nuc_GPU(const, ion, params, beta, field, eSpin, nSpin);
            end
            % end

            if accpRate(time, tt, ff) <= 0.0001 ||...
                    abs( dE(time)/sum(E_si(:,tt,ff) + E_int(:,tt,ff), 1) ) <= params.convg
                trap_check = trap_check + 1;
                if trap_check >= 20
                    break
                end
            end
            Et_temp(time) = sum(E_si(:,tt,ff) + E_int(:,tt,ff)); % record global energy
            if mod(time, 200) == 0; fprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                    accpRate(time, tt, ff)*100, time); end % checkpoint
        end
        spin_e(:,:,tt,ff) = eSpin;
        spin_n(:,:,tt,ff) = nSpin;
        E_tot(:, tt, ff) = Et_temp;
        dEt(:, tt, ff, 1) = dE;
        fprintf('Thermalization complete, total iteration: %.3e\n', time);
    end
end
% Retrieve results from GPU
init.eSpinT = gather(spin_e); % intermediate electronic spin configuration after thermalization
init.nSpinT = gather(spin_n); % intermediate nuclear spin configuration after thermalization
init.Etot = gather(E_tot); % global energy change during thermalization
init.dEt = gather(dEt); % history of energy change during thermalization
init.aRate = gather(accpRate); % history of acceptance rate during thermalization
init.Esi = E_si; % single-ion energy
init.Eint = E_int;
end