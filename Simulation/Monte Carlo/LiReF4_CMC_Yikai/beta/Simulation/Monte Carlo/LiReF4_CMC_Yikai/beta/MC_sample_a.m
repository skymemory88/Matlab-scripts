function [Mx, My, Mz, coef, eSpin, nSpin, E_si, E_int] = MC_sample_a(const, ion, params, E_si, hamI, basis, coef, eSpin, nSpin)
persistent time
if isempty(time); time = 0; end

if params.useGPU && ~isempty(params.gpuState)
    gpuState = params.gpuState;
else
    gpuState = [];
end
Mx = double.empty(length(params.temp), size(params.field, 2),0); % <Jx>
My = double.empty(length(params.temp), size(params.field, 2),0); % <Jy>
Mz = double.empty(length(params.temp), size(params.field, 2),0); % <Jz>
E_int = double.empty(length(params.temp), size(params.field,2),0); % interaction energy

% Initialize autocorrelation tracking parameters
autocorr_window = 100; % Window size for autocorrelation calculation
autocorr_threshold = 5; % If tau > threshold, use cluster updates (lowered from 10)
cluster_interval_min = 50; % Minimum steps between cluster attempts
cluster_interval_max = 500; % Maximum steps between cluster attempts
acceptance_threshold = 0.1; % Use clusters if acceptance < 10% (raised from 1%)
force_cluster_every = 1000; % Force a cluster attempt every N steps regardless

if length(params.temp) <= size(params.field,2)
    temp = params.temp;
    beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
    beta(beta==inf) = 1e9; % soften the divergence

    for tt = 1:length(temp)
        parfor ff = 1:size(params.field,2)
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration

            % Local variables for autocorrelation tracking (parfor compatible)
            local_history = zeros(1, autocorr_window);
            history_idx = 0;
            local_tau = 1;
            last_cluster_step = 0;
            cluster_interval = cluster_interval_min;

            for tIntv = 1:params.mIntv
                % Regular single-spin updates
                [~, accpRate, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize_a(ion, const, params,...
                    beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);

                % Thermalize nuclear spins
                if params.hyp
                    nSpinT = therm_nuc_a(const, beta(tt), ion, params, ff, eSpinT, nSpinT);
                end

                % Track observable for autocorrelation (using magnetization magnitude)
                current_mag = sqrt(mean(eSpinT(:,1))^2 + mean(eSpinT(:,2))^2 + mean(eSpinT(:,3))^2);

                % Update circular buffer
                history_idx = mod(history_idx, autocorr_window) + 1;
                local_history(history_idx) = current_mag;

                % Calculate autocorrelation time periodically after buffer is full
                if mod(tIntv, 50) == 0 && tIntv >= autocorr_window
                    % Reorder circular buffer for autocorrelation calculation
                    if history_idx == autocorr_window
                        ordered_history = local_history;
                    else
                        ordered_history = [local_history(history_idx+1:end), local_history(1:history_idx)];
                    end

                    local_tau = computeAutocorrTime(ordered_history);

                    % Adaptive cluster interval based on autocorrelation
                    if local_tau > autocorr_threshold
                        % High autocorrelation - more frequent clusters
                        cluster_interval = max(cluster_interval_min, ...
                            round(cluster_interval * 0.9));
                    else
                        % Low autocorrelation - less frequent clusters
                        cluster_interval = min(cluster_interval_max, ...
                            round(cluster_interval * 1.1));
                    end
                end

                % Decide whether to attempt cluster update
                steps_since_cluster = tIntv - last_cluster_step;
                should_cluster = (local_tau > autocorr_threshold && ...
                    steps_since_cluster >= cluster_interval) || ...
                    (accpRate < acceptance_threshold && steps_since_cluster >= cluster_interval_min) || ...
                    (steps_since_cluster >= force_cluster_every); % Force periodic attempts

                % Diagnostic output on first few decisions
                if tIntv <= 200 && mod(tIntv, 50) == 0 && isfield(params, 'verbose') && params.verbose
                    fprintf('Step %d: τ=%.2f (threshold=%.1f), accpRate=%.3f (threshold=%.2f), cluster in %d steps\n', ...
                        tIntv, local_tau, autocorr_threshold, accpRate, acceptance_threshold, ...
                        cluster_interval - steps_since_cluster);
                end

                if should_cluster
                    % Attempt cluster update
                    [cSize, de, cAccept, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = ...
                        thermalize_cluster_a(ion, const, params, beta(tt), ...
                        E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), ...
                        basis(:,:,tt,ff), eSpinT, nSpinT, tt, ff);

                    last_cluster_step = tIntv;

                    % Always log cluster attempts for debugging
                    if isfield(params, 'verbose') && params.verbose
                        fprintf('T=%.3f: Step %d, Cluster attempt - size=%d, ΔE=%.4f, %s (τ=%.1f, accpRate=%.3f)\n', ...
                            1/(const.kB*beta(tt)*const.J2meV), tIntv, cSize, de, ...
                            sprintf('%s', cAccept*'accepted' + ~cAccept*'rejected'), ...
                            local_tau, accpRate);
                    end
                end

                % Periodic progress report
                if isfield(params, 'verbose') && params.verbose && mod(tIntv, 200) == 0
                    fprintf('T=%.3f, B=%.3f: Step %d, τ=%.1f, cluster interval=%d\n', ...
                        temp(tt), vecnorm(params.field(:,ff)), tIntv, local_tau, cluster_interval);
                end
            end

            % Record measurements
            Mx(tt,ff,1) = mean(eSpinT(:,1)); % <Jx>
            My(tt,ff,1) = mean(eSpinT(:,2)); % <Jy>
            Mz(tt,ff,1) = mean(eSpinT(:,3)); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;

            % Calculate interaction energy (GPU or CPU)
            if params.useGPU && ~isempty(params.gpuState)
                gpuState.updateSpins(eSpinT, tt, ff);
                E_int(tt,ff,1) = dipSum_gpu_a(gpuState, tt, ff);
            else
                E_int(tt,ff,1) = dipSum_a(params, const.gfac, eSpinT);
            end
        end
    end
    time = time + params.mIntv;
else
    % Field sweep version
    for ff = 1:size(params.field,2)
        temp = params.temp;
        beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
        beta(beta==inf) = 1e9; % soften the divergence

        parfor tt = 1:length(temp)
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration

            % Local tracking for autocorrelation
            local_history = zeros(1, autocorr_window);
            history_idx = 0;
            local_tau = 1;
            last_cluster_step = 0;
            cluster_interval = cluster_interval_min;

            for tIntv = 1:params.mIntv
                % Regular updates
                [~, accpRate, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize_a(ion, const, params,...
                    beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);

                if params.hyp
                    nSpinT = therm_nuc_a(const, beta(tt), ion, params, ff, eSpinT, nSpinT);
                end

                % Track magnetization for autocorrelation
                current_mag = sqrt(mean(eSpinT(:,1))^2 + mean(eSpinT(:,2))^2 + mean(eSpinT(:,3))^2);

                % Update circular buffer
                history_idx = mod(history_idx, autocorr_window) + 1;
                local_history(history_idx) = current_mag;

                % Calculate autocorrelation
                if mod(tIntv, 50) == 0 && tIntv >= autocorr_window
                    % Reorder circular buffer
                    if history_idx == autocorr_window
                        ordered_history = local_history;
                    else
                        ordered_history = [local_history(history_idx+1:end), local_history(1:history_idx)];
                    end

                    local_tau = computeAutocorrTime(ordered_history);

                    if local_tau > autocorr_threshold
                        cluster_interval = max(cluster_interval_min, round(cluster_interval * 0.9));
                    else
                        cluster_interval = min(cluster_interval_max, round(cluster_interval * 1.1));
                    end
                end

                % Cluster update decision
                steps_since_cluster = tIntv - last_cluster_step;
                should_cluster = (local_tau > autocorr_threshold && ...
                    steps_since_cluster >= cluster_interval) || ...
                    (accpRate < 0.01 && steps_since_cluster >= cluster_interval_min);

                if should_cluster
                    [~, ~, ~, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = ...
                        thermalize_cluster_a(ion, const, params, beta(tt), ...
                        E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), ...
                        basis(:,:,tt,ff), eSpinT, nSpinT, tt, ff);

                    last_cluster_step = tIntv;
                end
            end

            % Record measurements
            Mx(tt,ff,1) = mean(eSpinT(:,1)); % <Jx>
            My(tt,ff,1) = mean(eSpinT(:,2)); % <Jy>
            Mz(tt,ff,1) = mean(eSpinT(:,3)); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;

            % Calculate interaction energy (GPU or CPU)
            if params.useGPU && ~isempty(params.gpuState)
                gpuState.updateSpins(eSpinT, tt, ff);
                E_int(tt,ff,1) = dipSum_gpu_a(gpuState, tt, ff);
            else
                E_int(tt,ff,1) = dipSum_a(params, const.gfac, eSpinT);
            end
        end
    end
    time = time + params.mIntv;
end
params.gpuState = gpuState;
end

% Helper function to compute integrated autocorrelation time
function tau = computeAutocorrTime(data)
% Compute integrated autocorrelation time using automated windowing
% Based on Sokal's adaptive truncation method

n = length(data);
data = data - mean(data); % Center the data

% Compute autocorrelation function
c0 = sum(data.^2) / n;
if c0 == 0
    tau = 1;
    return;
end

tau_int = 0.5; % Start with C(0)/2

% Adaptive windowing following Sokal
for t = 1:min(n-1, n/4)
    ct = sum(data(1:n-t) .* data(t+1:n)) / (n-t) / c0;
    tau_int = tau_int + ct;

    % Sokal's criterion: stop when t >= c * tau_int
    if t >= 6 * tau_int
        break;
    end
end

tau = max(1, 2*tau_int); % Factor of 2 for full autocorrelation time
end