function [spin_e, spin_n, coef, Esi, Edip, E_nuc, dEt, accpRate] = equilibrate_a(const, ion, params, Esi, Edip, hamI, coef, basis, spin_e, E_nuc, spin_n)
t_lim = params.tEQ; % equilibration time limit
accpRate = zeros(t_lim, length(params.temp), size(params.field, 2)); % acceptance rate tracker
dEt = double.empty(t_lim, length(params.temp), size(params.field, 2), 0); % global energy update trace
for tt = 1:length(params.temp)
    parfor ff = 1:size(params.field,2)
        worker = getCurrentTask; % obtain CPU core ID  
        if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
            cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
        else
            cID = '';
        end
        eSpin = squeeze(spin_e(:,:,tt,ff)); % initial electronic spin configuration
        nSpin = squeeze(spin_n(:,:,tt,ff)); % initial nuclear spin configuration
        Enn = squeeze(E_nuc(:,tt,ff)); % initial nuclear spin energy
        fprintf([cID sprintf('Temperature: %.2f, Field: %.2f.\n', params.temp(tt), vecnorm(params.field(:,ff)))])

        trap_check = 0; % critical slowing-down token
        convgTkr = 0; % convergence tracker
        dE = zeros(t_lim,1);
        for time = 1:t_lim         
            % thermalize the electronic spins
            [dE(time), accpRate(time,tt,ff), Esi(:,tt,ff), Edip(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_a(ion, const, params,...
                tt, Esi(:,tt,ff), Edip(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpin, nSpin);
            if accpRate(time, tt, ff) <= 0.05 || abs( dE(time)/sum(Esi(:,tt,ff) + Edip(:,tt,ff)) ) <= params.convg
                trap_check = trap_check + 1;
                if trap_check >= 5 % after 5 consecutive trapped random walk, invoke cluster update
                    [cSize, de, aRate, Esi(:,tt,ff), Edip(:,tt,ff), coef(:,:,tt,ff), eSpin] = thermalize_cluster_a(ion, const, params,...
                        tt, Esi(:,tt,ff), Edip(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpin, nSpin);
                    dE(time) = de * aRate;
                    fprintf([cID sprintf('Cluster size: %u, acceptance: %u, iteration: %.3e\n', cSize, aRate, time)]);
                    trap_check = 0; % reset the tracker
                    if aRate == 0; convgTkr = convgTkr + 1; end
                    if convgTkr >= 5; break; end
                    % break; % when cluster update is not used
                end
            end          
            % thermalize the nuclear spins
            if params.hyp
                [nSpin, Enn] = therm_nuc_a(const, ion, params, tt, ff, eSpin, nSpin, Enn);
            end
            if mod(time, 200) == 0; fprintf([cID sprintf('Current acceptance rate: %.2f, iteration: %.3e\n', accpRate(time, tt, ff), time)]); end % checkpoint
        end
        spin_e(:,:,tt,ff) = eSpin;
        spin_n(:,:,tt,ff) = nSpin;
        dEt(:, tt, ff, 1) = dE;
        E_nuc(:,tt,ff) = Enn;
        fprintf([cID sprintf('Thermalization complete, total iteration: %.3e\n', time)]);
    end
end
end