function [init, coef, E_si] = equilibrate_GPU(const, ion, params, E_si, hamI, coef, basis, eSpinT, nSpinT)
% Convert parameters and initial arrays to GPU arrays
t_lim = params.tEQ; % equilibration time limit
accpRate = zeros(params.tEQ, length(params.temp), size(params.field, 2)); % acceptance rate tracker
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
% Main equilibration loop
for tt = 1:length(params.temp)
    beta = gpuArray(1 / (const.J2meV * const.kB * params.temp(tt)));
    for ff = 1:size(params.field,2)
        field = gpuArray(params.field(:,ff));
        eSpin = squeeze(eSpinT(:,:,tt,ff)); % initial electronic spin configuration
        nSpin = squeeze(nSpinT(:,:,tt,ff)); % initial nuclear spin configuration
        fprintf('Temperature: %.3f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))

        trap_check = 0; % critical slowing-down token
        % convgTkr = 0; % cluster update convergence tracker
        Et_si = zeros(t_lim,1);
        dE = zeros(t_lim,1);
        for time = 1:t_lim
            % thermalize the electronic spins
            [E_si(:,tt,ff), dE(time), coef(:,:,tt,ff), eSpin, accpRate(time,tt,ff)] = thermalize_GPU(const, ion,...
                params, E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin, beta);
            % thermalize the nuclear spins
            % if params.hyp && mod(time,20) == 0 % wot slow down nuclear spin update
            if params.hyp
                nSpin = therm_nuc_GPU(const, ion, params, beta, field, eSpin, nSpin);
            end

            Et_si(time) = sum(E_si(:,tt,ff)); % update global energy
            
            if mod(time, 200) == 0; fprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                    accpRate(time, tt, ff)*100, time); end % checkpoint
            if accpRate(time, tt, ff) <= 0.0001 ||...
                    abs( dE(time)/sum(E_si(:,tt,ff), 1) ) <= params.convg
                trap_check = trap_check + 1;
                if trap_check >= 20
                    disp('Convergence reached!')
                    break
                end
            end
            if time == t_lim; disp('Thermalization step limit reached!'); end
        end
        eSpinT(:,:,tt,ff) = eSpin;
        nSpinT(:,:,tt,ff) = nSpin;
        dEt(:, tt, ff, 1) = dE;
        fprintf('Thermalization complete, total iteration: %.3e\n', time);
    end
    % Annealing: use the final state of (tt+1) temperature as the starting configuration of (tt) temperature
    if tt < length(params.temp)
        eSpinT(:,:,tt+1,:) = eSpin;
        nSpinT(:,:,tt+1,:) = nSpin;
        E_si(:,tt+1,:) = E_si(:,tt,:);
    end
end
% Retrieve results from GPU
init.eSpinT = gather(eSpinT); % intermediate electronic spin configuration after thermalization
init.nSpinT = gather(nSpinT); % intermediate nuclear spin configuration after thermalization
init.dEt = gather(dEt); % history of energy change during thermalization
init.aRate = gather(accpRate); % history of acceptance rate during thermalization
init.Esi = E_si; % single-ion energy
end