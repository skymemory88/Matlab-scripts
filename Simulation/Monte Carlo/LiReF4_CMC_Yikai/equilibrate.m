function [spin_e, spin_n, coef, E_si, dEt, accpRate] = equilibrate(const, ion, params, hamI, coef, basis, E_si, spin_e, spin_n)
t_lim = params.tEQ; % equilibration time limit
pt_intv = params.pt_intv; % parallel tempering intervals
pt_steps = int32(t_lim / pt_intv);
if length(params.temp) <= size(params.field,2)
    temp = params.temp;
    accpRate = zeros(t_lim, length(temp), size(params.field, 2)); % acceptance rate tracker
    dEt = double.empty(t_lim, length(temp), size(params.field, 2), 0); % global energy update trace
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
                [dE(step), accpRate(step,tt,ff), E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize(ion, const, params,...
                    bet, E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);

                % thermalize the nuclear spins
                % if params.hyp && mod(time,20) == 0 % slow down nuclear spin update
                if params.hyp
                    [nSpin] = therm_nuc(const, bet, ion, params, ff, eSpin, nSpin);
                end

                if accpRate(step, tt, ff) <= 0.0001 || abs( dE(step)/sum(E_si(:,tt,ff), 1) ) <= params.convg
                    trap_check = trap_check + 1;
                    % if trap_check >= 5 && trap_check < 10 % after 5 consecutive trapped random walk, invoke cluster update
                    %     [cSize, de, aRate, Esi(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_cluster(ion, const, params,...
                    %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpin, nSpin);
                    %     dE(time) = de * aRate;
                    %     fprintf([cID sprintf('Cluster size: %1$u, acceptance: %2$u, iteration: %3$.3e\n', cSize, aRate, time)]);
                    %     % if aRate == 0; convgTkr = convgTkr + 1; end
                    % end
                    if trap_check >= 20
                        break
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
    dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
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
                    [dE(ticks, step, tt), aRate(ticks,step,tt), E_si(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize(ion, const, params,...
                        beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);
                    % thermalize the nuclear spins
                    if params.hyp
                        [nSpin] = therm_nuc(const, beta(tt), ion, params, ff, eSpin, nSpin);
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
                end
                Et(tt,1) = sum(E_si(:,tt,ff),1) + dipSum(const.gfac, params.pos, eSpin);
                spin_e(:,:,tt,ff) = eSpin;
                spin_n(:,:,tt,ff) = nSpin;
            end

            % Parallel tempering algorithm
            if step < pt_steps % skip the last round
                [~, tIdx] = sort(temp); % indices of the temperatures in ascending order
                for tt = 1+rem(step,2):2:length(temp)-1 % iterate over neighbouring temperatures
                    Ediff = Et( tIdx(tt+1) ) - Et( tIdx(tt) );
                    crit = rand();
                    prob = min(1, Ediff * (beta( tIdx(tt+1) ) - beta( tIdx(tt) )));
                    if prob >= crit
                        [temp( tIdx(tt) ), temp( tIdx(tt+1) )] = deal(temp( tIdx(tt+1) ), temp( tIdx(tt) ));
                        fprintf('Swapped configuration between T = %.3f and T = %.3f.\n', temp(tIdx(tt)), temp(tIdx(tt+1)));
                    end
                end
                beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
                beta(beta==inf) = 1e9; % soften the divergence
            end
        end
        [~, tIdx] = sort(temp); % indices of the temperatures in ascending order
        spin_e(:,:,:,ff) = spin_e(:,:,tIdx,ff);
        spin_n(:,:,:,ff) = spin_n(:,:,tIdx,ff);
        dEt(:, :, ff, 1) = reshape(dE, pt_steps * pt_intv, length(temp));
        accpRate(:, :, ff, 1) = reshape(aRate, pt_steps * pt_intv, length(temp));
    end
end