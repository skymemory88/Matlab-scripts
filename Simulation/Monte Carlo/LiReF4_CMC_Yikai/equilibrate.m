function [spin, coef, en0, dEt, accpRate] = equilibrate(const, ion, params, en0, hamI, coef, basis, spin0)
pos = params.pos;
t_lim = params.tEQ; % equilibration time limit

accpRate = zeros(t_lim, length(params.temp), size(params.field, 2)); % acceptance rate tracker
spin = double.empty(size(pos,1), 3, length(params.temp), size(params.field, 2), 0); % final spin config
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
for tt = 1:length(params.temp)
    parfor ff = 1:size(params.field,2)
        worker = getCurrentTask; % obtain CPU core ID  
        if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
            cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
        else
            cID = '';
        end
       
        wav = coef(:,:,tt,ff); % wavfunction
        bas = basis(:,:,tt,ff); % Ising basis
        fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))])

        spin_i = squeeze(spin0(:,:,tt,ff)); % initial spin configuration
        trap_check = 0; % critical slowing-down token
        convgTkr = 0; % convergence tracker
        dE = zeros(t_lim,1);
        for time = 1:t_lim
            % thermalize the spin configuration
            [dE(time), accpRate(time,tt,ff), en0(:,tt,ff), wav, spin_i] = thermalize(ion, const, params, params.temp(tt), en0(:,tt,ff), hamI(:,:,tt,ff), wav, bas, spin_i);
            if accpRate(time, tt, ff) <= 0.05 || abs(dE(time)/sum(en0(:,tt,ff))) <= params.convg
                trap_check = trap_check + 1;
                if trap_check >= 5 % after 5 consecutive trapped random walk, invoke cluster update
                    [cSize, de, aRate, spin_i, wav, en0(:,tt,ff)] = thermalize_cluster(ion, const, params, params.temp(tt), en0(:,tt,ff), hamI(:,:,tt,ff), wav, spin_i);
                    dE(time) = de * aRate;
                    fprintf([cID sprintf('Cluster size: %u, acceptance: %u, iteration: %.3e\n', cSize, aRate, time)]);
                    trap_check = 0; % reset the tracker
                    if aRate == 0; convgTkr = convgTkr + 1; end
                    if convgTkr >= 10; break; end
                    % break;
                end
            end
            if mod(time, 200) == 0; fprintf([cID sprintf('Current acceptance rate: %.2f, iteration: %.3e\n', accpRate(time, tt, ff), time)]); end % checkpoint
        end
        dEt(:, tt, ff, 1) = dE;
        coef(:,:,tt,ff) = wav;
        spin(:,:,tt,ff,1) = spin_i;
        fprintf([cID sprintf('Thermalization complete, total iteration: %.3e\n', time)]);
    end
end
end