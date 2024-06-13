function [spin_e, spin_n, coef, E_si, E_int, E_tot, dEt, accpRate] = equilibrate_a(const, ion, params, hamI, coef, basis, E_si, E_int, E_tot, spin_e, spin_n)
t_lim = params.tEQ; % equilibration time limit
accpRate = zeros(t_lim, length(params.temp), size(params.field, 2)); % acceptance rate tracker
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
if length(params.temp) <= size(params.field,2)
    for tt = 1:length(params.temp)
        if params.temp(tt) > 0
            beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
        else
            beta = 1e9; % approximate infinity at zero temperature
        end
        parfor ff = 1:size(params.field,2)
            worker = getCurrentTask; % obtain CPU core ID
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
            nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
            fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))])

            trap_check = 0; % critical slowing-down token
            % convgTkr = 0; % cluster update convergence tracker
            Et_temp = zeros(t_lim,1);
            dE = zeros(t_lim,1);
            for time = 1:t_lim
                % thermalize the electronic spins
                [dE(time), accpRate(time,tt,ff), E_si(:,tt,ff), E_int(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_a(ion, const, params,...
                    beta, E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);

                % thermalize the nuclear spins
                % if params.hyp && mod(time,20) == 0 % slow down nuclear spin update
                if params.hyp
                    [nSpin] = therm_nuc_a(const, beta, ion, params, ff, eSpin, nSpin);
                end

                if accpRate(time, tt, ff) <= 0.0001 || abs( dE(time)/sum(E_si(:,tt,ff) + E_int(:,tt,ff), 1) ) <= params.convg
                    trap_check = trap_check + 1;
                    % if trap_check >= 5 && trap_check < 10 % after 5 consecutive trapped random walk, invoke cluster update
                    %     [cSize, de, aRate, Esi(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_cluster_a(ion, const, params,...
                    %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpin, nSpin);
                    %     dE(time) = de * aRate;
                    %     fprintf([cID sprintf('Cluster size: %1$u, acceptance: %2$u, iteration: %3$.3e\n', cSize, aRate, time)]);
                    %     % if aRate == 0; convgTkr = convgTkr + 1; end
                    % end
                    if trap_check >= 20
                        break
                    end
                end
                Et_temp(time) = sum(E_si(:,tt,ff) + E_int(:,tt,ff)); % record global energy
                if mod(time, 200) == 0; disp( strcat(cID, sprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                        accpRate(time, tt, ff)*100, time)) ); end % checkpoint
            end
            spin_e(:,:,tt,ff) = eSpin;
            spin_n(:,:,tt,ff) = nSpin;
            E_tot(:, tt, ff) = Et_temp;
            dEt(:, tt, ff, 1) = dE;
            fprintf([cID sprintf('Thermalization complete, total iteration: %.3e\n', time)]);
        end
    end
else
    for ff = 1:size(params.field,2)
        parfor tt = 1:length(params.temp)
        % for tt = 1:length(params.temp) % for debugging
            if params.temp(tt) > 0
                beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
            else
                beta = 1e9; % approximate infinity at zero temperature
            end
            worker = getCurrentTask; % obtain CPU core ID
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
            nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
            fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))])

            trap_check = 0; % critical slowing-down token
            % convgTkr = 0; % cluster update convergence tracker
            Et_temp = zeros(t_lim,1);
            dE = zeros(t_lim,1);
            for time = 1:t_lim
                % thermalize the electronic spins
                [dE(time), accpRate(time,tt,ff), E_si(:,tt,ff), E_int(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_a(ion, const, params,...
                    beta, E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);

                % thermalize the nuclear spins
                % if params.hyp && mod(time,20) == 0 % slow down nuclear spin update
                if params.hyp
                    [nSpin] = therm_nuc_a(const, beta, ion, params, ff, eSpin, nSpin);
                end

                if accpRate(time, tt, ff) <= 0.0001 || abs( dE(time)/sum(E_si(:,tt,ff)) ) <= params.convg
                    trap_check = trap_check + 1;
                    % if trap_check >= 5 && trap_check < 10 % after 5 consecutive trapped random walk, invoke cluster update
                    %     [cSize, de, aRate, Esi(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_cluster_a(ion, const, params,...
                    %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpin, nSpin);
                    %     dE(time) = de * aRate;
                    %     fprintf([cID sprintf('Cluster size: %1$u, acceptance: %2$u, iteration: %3$.3e\n', cSize, aRate, time)]);
                    %     % if aRate == 0; convgTkr = convgTkr + 1; end
                    % end
                    if trap_check >= 20
                        break
                    end
                end
                Et_temp(time) = sum(E_si(:,tt,ff) + E_int(:,tt,ff)); % record global energy
                if mod(time, 200) == 0; disp( strcat(cID, sprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                        accpRate(time, tt, ff)*100, time)) ); end % checkpoint
            end
            spin_e(:,:,tt,ff) = eSpin;
            spin_n(:,:,tt,ff) = nSpin;
            E_tot(:, tt, ff) = Et_temp;
            dEt(:, tt, ff, 1) = dE;
            fprintf([cID sprintf('Thermalization complete, total iteration: %.3e\n', time)]);
        end
    end
end
end