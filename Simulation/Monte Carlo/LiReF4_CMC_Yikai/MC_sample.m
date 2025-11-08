function [Mx, My, Mz, coef, eSpin, nSpin, E_si, E_int] = MC_sample(const, ion, params, E_si, hamI, basis, coef, eSpin, nSpin)
persistent time
if isempty(time); time = 0;end
ticker = time;
Mx = double.empty(length(params.temp), size(params.field, 2),0); % <Jx>
My = double.empty(length(params.temp), size(params.field, 2),0); % <Jy>
Mz = double.empty(length(params.temp), size(params.field, 2),0); % <Jz>
E_int = double.empty(length(params.temp), size(params.field,2),0); % interaction energy
if length(params.temp) <= size(params.field,2)
    temp = params.temp;
    beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
    beta(beta==inf) = 1e9; % soften the divergence
    for tt = 1:length(temp)
        parfor ff = 1:size(params.field,2)
        % for ff = 1:size(params.field,2)
            worker = getCurrentTask; % For parallel execution only
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration
            for tIntv = 1:params.mIntv
                [~, ~, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize(ion, const, params,...
                    beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);
                % thermalize the nuclear spins
                if params.hyp
                    nSpinT = therm_nuc(const, beta(tt), ion, params, ff, eSpinT, nSpinT);
                end
                % checkpoint
                if mod(tIntv, ceil(params.mIntv/2)) == 0
                    fprintf([cID sprintf('Current iteration: %.3e \n', tIntv + ticker)]);
                end
            end
            % % try cluster update for five times after every data acquisition
            % for trial = 1:5
            %     [~, ~, ~, Esi(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize_cluster(ion, const, params,...
            %         tt, Esi(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), eSpinT, nSpinT);
            % end

            % Record measurements
            Mx(tt,ff,1) = rms(eSpinT(:,1),1); % <Jx>
            My(tt,ff,1) = rms(eSpinT(:,2),1); % <Jy>
            Mz(tt,ff,1) = rms(eSpinT(:,3),1); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;
            E_int(tt,ff,1) = dipSum(const.gfac, params.pos, eSpin(:,:,tt,ff));
        end
    end
    time = time + params.mIntv;
else
    for ff = 1:size(params.field,2)
    temp = params.temp;
    beta = 1 ./ (const.J2meV * const.kB * temp); % 1/kBT [meV]
    beta(beta==inf) = 1e9; % soften the divergence
        parfor tt = 1:length(temp)
        % for tt = 1:length(temp)
            worker = getCurrentTask; % For parallel execution only
            if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
                cID = sprintf('Core %d: ', worker.ID); % For parallel execution only
            else
                cID = '';
            end
            eSpinT = squeeze(eSpin(:,:,tt,ff)); % electronic spin configuration
            nSpinT = squeeze(nSpin(:,:,tt,ff)); % nuclear spin configuration

            for tIntv = 1:params.mIntv
                [~, ~, E_si(:,tt,ff), coef(:,:,tt,ff), eSpinT] = thermalize(ion, const, params,...
                    beta(tt), E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);
                % thermalize the nuclear spins
                if params.hyp
                    nSpinT = therm_nuc(const, beta(tt), ion, params, ff, eSpinT, nSpinT);
                end
                % checkpoint
                if mod(tIntv, ceil(params.mIntv/2)) == 0
                    fprintf([cID sprintf('Current iteration: %.3e\n', tIntv)]);
                end
            end

            % Record measurements
            Mx(tt,ff,1) = rms(eSpinT(:,1),1); % <Jx>
            My(tt,ff,1) = rms(eSpinT(:,2),1); % <Jy>
            Mz(tt,ff,1) = rms(eSpinT(:,3),1); % <Jz>
            eSpin(:,:,tt,ff) = eSpinT;
            nSpin(:,:,tt,ff) = nSpinT;
            E_int(tt,ff,1) = dipSum(const.gfac, params.pos, eSpin(:,:,tt,ff));
        end
    end
    time = time + params.mIntv;
end
end