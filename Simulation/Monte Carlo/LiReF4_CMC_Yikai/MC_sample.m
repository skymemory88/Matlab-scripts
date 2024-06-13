function [Mx, My, Mz, coef, eSpin, nSpin, E_si, E_int] = MC_sample(const, ion, params, E_si, E_int, hamI, basis, coef, eSpin, nSpin)
persistent time
if isempty(time); time = 0;end
ticker = time;
Mx = double.empty(1, length(params.temp), size(params.field, 2),0); % <Jx>
My = double.empty(1, length(params.temp), size(params.field, 2),0); % <Jy>
Mz = double.empty(1, length(params.temp), size(params.field, 2),0); % <Jz>
if length(params.temp) <= size(params.field,2)
    for tt = 1:length(params.temp)
        if params.temp(tt) > 0
            beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
        else
            beta = 1e9; % approximate infinity at zero temperature
        end
        parfor ff = 1:size(params.field,2)
            worker = getCurrentTask; % For parallel execution only
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration
            for tInt = 1:params.mIntv
                [~, ~, E_si(:,tt,ff), E_int(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize(ion, const, params,...
                    beta, E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);
                % thermalize the nuclear spins
                if params.hyp
                    nSpinT = therm_nuc(const, beta, ion, params, ff, eSpinT, nSpinT);
                end
                % checkpoint
                if mod(tInt, ceil(params.mIntv/2)) == 0
                    fprintf([cID sprintf('Current iteration: %.3e \n', tInt + ticker)]);
                end
            end
            % % try cluster update for five times after every data acquisition
            % for trial = 1:5
            %     [~, ~, ~, Esi(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize_cluster(ion, const, params,...
            %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpinT, nSpinT);
            % end

            % Measurement
            Mx(1, tt,ff,1) = rms(eSpinT(:,1),1); % <Jx>
            My(1, tt,ff,1) = rms(eSpinT(:,2),1); % <Jy>
            Mz(1, tt,ff,1) = rms(eSpinT(:,3),1); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;
        end
    end
    time = time + params.mIntv;
else
    for ff = 1:size(params.field,2)
        parfor tt = 1:length(params.temp)
            if params.temp(tt) > 0
                beta = 1/ (const.J2meV * const.kB * params.temp(tt)); % [meV]
            else
                beta = 1e9; % approximate infinity at zero temperature
            end
            worker = getCurrentTask; % For parallel execution only
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration

            for tInt = 1:params.mIntv
                [~, ~, E_si(:,tt,ff), E_int(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize(ion, const, params,...
                    beta, E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);
                % thermalize the nuclear spins
                if params.hyp
                    nSpinT = therm_nuc(const, beta, ion, params, ff, eSpinT, nSpinT);
                end
                % checkpoint
                if mod(tInt, ceil(params.mIntv/2)) == 0
                    fprintf([cID sprintf('Current iteration: %.3e\n', tInt)]);
                end
            end
            % % try cluster update for five times after every data acquisition
            % for trial = 1:5
            %     [~, ~, ~, Esi(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize_cluster(ion, const, params,...
            %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpinT, nSpinT);
            % end

            % Measurement
            Mx(1, tt,ff,1) = rms(eSpinT(:,1),1); % <Jx>
            My(1, tt,ff,1) = rms(eSpinT(:,2),1); % <Jy>
            Mz(1, tt,ff,1) = rms(eSpinT(:,3),1); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;
        end
    end
    time = time + params.mIntv;
end
end