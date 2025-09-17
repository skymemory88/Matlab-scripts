function [spin_e, spin_n, coef, E_si, dEt, accpRate] = equilibrate_a(const, ion, params, hamI, coef, basis, E_si, spin_e, spin_n)
t_lim = params.tEQ; % equilibration time limit
pt_intv = params.pt_intv; % parallel tempering intervals
pt_steps = int32(t_lim / pt_intv);
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
if length(params.temp) <= size(params.field,2)
    temp = params.temp;
    accpRate = zeros(t_lim, length(temp), size(params.field, 2)); % acceptance rate tracker
    beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
    beta(beta==inf) = 1e9; % soften the divergence
    for tt = 1:length(temp)
        bet = beta(tt);
        parfor ff = 1:size(params.field,2)
            % for ff = 1:size(params.field,2)
            worker = getCurrentTask; % obtain CPU core ID
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
            nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
            fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', temp(tt), vecnorm(params.field(:,ff)))])

            trap_check = 0; % critical slowing-down token
            % convgTkr = 0; % cluster update convergence tracker
            dE = zeros(t_lim,1);
            for step = 1:t_lim
                % thermalize the electronic spins
                [dE(step), accpRate(step,tt,ff), E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_a(ion, const, params,...
                    bet, E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);

                % thermalize the nuclear spins
                % if params.hyp && mod(time,20) == 0 % slow down nuclear spin update
                if params.hyp
                    [nSpin] = therm_nuc_a(const, bet, ion, params, ff, eSpin, nSpin);
                end

                if trap_check >= 5 && trap_check < 10 % after 5 consecutive trapped random walks
                    % Attempt cluster update
                    [cSize, de, aRate, E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = ...
                        thermalize_cluster_a(ion, const, params, bet, E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin, tt, ff);

                    dE(step) = de * aRate; % Record energy change
                    if aRate > 0
                        fprintf([cID sprintf('Cluster update: size=%d, ΔE=%.4f meV, accepted, iteration: %d\n', ...
                            cSize, de, step)]);
                        trap_check = 0; % Reset trap counter on successful cluster update
                    else
                        fprintf([cID sprintf('Cluster update: size=%d, rejected, iteration: %d\n', ...
                            cSize, step)]);
                    end
                end
                if mod(step, 200) == 0; disp( strcat(cID, sprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                        accpRate(step, tt, ff)*100, step)) ); end % checkpoint
            end
            spin_e(:,:,tt,ff) = eSpin;
            spin_n(:,:,tt,ff) = nSpin;
            dEt(:, tt, ff, 1) = dE;
            fprintf([cID sprintf('Thermalization complete, total iteration: %.3e\n', step)]);
        end
    end
else
    gpuState = params.gpuState;
    accpRate = zeros(t_lim, length(params.temp), size(params.field, 2)); % acceptance rate tracker
    for ff = 1:size(params.field,2)
        temp = params.temp;
        beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
        beta(beta==inf) = 1e9; % soften the divergence
        dE = zeros(pt_intv, pt_steps, length(temp)); % change in global energy
        aRate = zeros(pt_intv, pt_steps, length(temp)); % acceptance rate
        Et = double.empty(length(temp), 0); % global energy update trace
        for step = 1:pt_steps
            parfor tt = 1:length(temp)
                % for tt = 1:length(temp) % for debugging
                worker = getCurrentTask; % obtain CPU core ID
                if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                    cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
                else
                    cID = '';
                end
                eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
                nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
                trap_check = 0;
                for ticks = 1:pt_intv
                    % thermalize the electronic spins
                    [dE(ticks, step, tt), aRate(ticks,step,tt), E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_a(ion, const, params,...
                        beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);

                    % thermalize the nuclear spins
                    if params.hyp
                        [nSpin] = therm_nuc_a(const, beta(tt), ion, params, ff, eSpin, nSpin);
                    end
                    % checkpoint
                    if mod(ticks, floor(pt_intv/3)) == 0
                        fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', temp(tt), vecnorm(params.field(:,ff)))])
                        disp( strcat(cID, sprintf('Current acceptance rate: %.2f%%, iteration: %.3e\n',...
                            aRate(ticks, step, tt)*100, ticks*(step + 1))) );
                    end
                    if aRate(ticks,step,tt) <= 0.0001 || abs( dE(ticks, step, tt)/sum(E_si(:,tt,ff), 1) ) <= params.convg
                        trap_check = trap_check + 1;
                        if trap_check >= params.trap_lim; break; end
                    end
                    
                    % Add cluster update for field-sweep scenario
                    if trap_check >= 5 && trap_check < 10 % after 5 consecutive trapped random walks
                        % Attempt cluster update
                        [cSize, de, aRate_cluster, E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = ...\
                            thermalize_cluster_a(ion, const, params, beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin, tt, ff);
                        
                        if aRate_cluster > 0
                            fprintf([cID sprintf('Field-sweep cluster update: size=%d, ΔE=%.4f meV, accepted, step=%d, tick=%d\\n', ...\
                                cSize, de, step, ticks)]);
                            trap_check = 0; % Reset trap counter on successful cluster update
                        else
                            fprintf([cID sprintf('Field-sweep cluster update: size=%d, rejected, step=%d, tick=%d\\n', ...\
                                cSize, step, ticks)]);
                        end
                    end
                end
                % Et(tt,1) = sum(E_si(:,tt,ff),1)...
                %     + dipSum_a(const.gfac, params.pos, eSpin, params.dpRng, params.nList, params.nVecs, params.nDists);
                spin_e(:,:,tt,ff) = eSpin;
                spin_n(:,:,tt,ff) = nSpin;
            end

            for tt = 1:length(temp)
                if params.useGPU && ~isempty(params.gpuState)
                    params.gpuState.updateSpins(squeeze(spin_e(:,:,tt,ff)), tt, ff);
                    E_dip = dipSum_gpu_a(params.gpuState, tt, ff);
                else
                    E_dip = dipSum_a(params, const.gfac, squeeze(spin_e(:,:,tt,ff)));
                end
                Et(tt) = sum(E_si(:,tt,ff)) + E_dip;
            end

            % Parallel tempering algorithm
            if step < pt_steps % skip the last round
                [~, tIdx] = sort(temp); % indices of the temperatures in ascending order
                for tt = 1+rem(step,2):2:length(temp)-1 % iterate over neighbouring temperatures
                    dBeta = beta(tIdx(tt)) - beta(tIdx(tt+1));
                    Ediff = Et(tIdx(tt+1)) - Et(tIdx(tt));
                    prob = min(1, exp(dBeta * Ediff));
                    crit = rand();
                    % if prob >= crit
                    %     [temp( tIdx(tt) ), temp( tIdx(tt+1) )] = deal(temp( tIdx(tt+1) ), temp( tIdx(tt) ));
                    %     fprintf('Swapped configuration between T = %.3f and T = %.3f.\n', temp(tIdx(tt)), temp(tIdx(tt+1)));
                    % end
                    if prob >= crit
                        [temp(tIdx(tt)), temp(tIdx(tt+1))] = deal(temp(tIdx(tt+1)), temp(tIdx(tt)));
                        fprintf('Swapped configuration between T = %.3f and T = %.3f.\n', temp(tIdx(tt)), temp(tIdx(tt+1)));
                    end
                end
                beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
                beta(beta==inf) = 1e9; % soften the divergence
            end
        end
        % sort the configuration in ascending temperatures
        [~, tIdx] = sort(temp); % indices of the temperatures in ascending order
        spin_e(:,:,:,ff) = spin_e(:,:,tIdx,ff);
        spin_n(:,:,:,ff) = spin_n(:,:,tIdx,ff);
        E_si(:,:,ff) = E_si(:,tIdx,ff); % added on 22.06.2025
        coef(:,:,:,ff) = coef(:,:,tIdx,ff); % added on 22.06.2025
        
        % update GPU spins
        if params.useGPU && ~isempty(params.gpuState)
            for tt = 1:length(params.temp)
                params.gpuState.updateSpins(squeeze(spin_e(:,:,tt,ff)), tt, ff);
            end
        end

        % record the MC update history
        dEt(:, :, ff, 1) = reshape(dE, pt_steps * pt_intv, length(temp));
        accpRate(:, :, ff) = reshape(aRate, pt_steps * pt_intv, length(temp));
    end
    params.gpuState = gpuState;
end