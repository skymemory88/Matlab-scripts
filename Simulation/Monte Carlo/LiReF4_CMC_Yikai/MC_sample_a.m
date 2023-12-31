function [Mx, My, Mz, wv, spin_e, spin_n, Esi, Edip, Enuc] = MC_sample_a(const, ion, params, Esi, Edip, hamI, basis, wv, spin_e, spin_n, Enuc)

Mx = double.empty(params.sIntv,length(params.temp), size(params.field, 2),0); % <Jx>
My = double.empty(params.sIntv,length(params.temp), size(params.field, 2),0); % <Jy>
Mz = double.empty(params.sIntv,length(params.temp), size(params.field, 2),0); % <Jz>
for tt = 1:length(params.temp)
    parfor ff = 1:size(params.field,2)
        worker = getCurrentTask; % For parallel execution only
        if parallel.internal.pool.isPoolThreadWorker || ~isempty(getCurrentJob)
            cID = fprintf('Core %d: ', worker.ID); % For parallel execution only
        else
            cID = '';
        end
        eSpinT = squeeze(spin_e(:,:,tt,ff)); % electronic spin configuration
        nSpinT = squeeze(spin_n(:,:,tt,ff)); % nuclear spin configuration
        Enn = squeeze(Enuc(:,tt,ff)); % nuclear spin energy

        time = 1; % calculation time tracker
        tSamp = 1; % sample size tracker
        mx = double.empty(params.sIntv,0);
        my = double.empty(params.sIntv,0);
        mz = double.empty(params.sIntv,0);
        while tSamp < params.sIntv
            for tInt = 1:params.mIntv
                [~, ~, Esi(:,tt,ff), Edip(:,tt,ff), wv(:,:,tt,ff), eSpinT] = thermalize_a(ion, const, params,...
                    tt, Esi(:,tt,ff), Edip(:,tt,ff), hamI(:,:,tt,ff), wv(:,:,tt,ff), basis(:,:,tt,ff), eSpinT, nSpinT);
                % thermalize the nuclear spins
                if params.hyp
                    [nSpinT, Enn] = therm_nuc_a(const, ion, params, tt, ff, eSpinT, nSpinT, Enn);
                end
                time = time + 1;
            end
            % Measurement
            mx(tSamp,1) = rms(eSpinT(:,1),1); % <Jx>
            my(tSamp,1) = rms(eSpinT(:,2),1); % <Jy>
            mz(tSamp,1) = rms(eSpinT(:,3),1); % <Jz>
            tSamp = tSamp + 1;
            if mod(time, 20) == 0
                fprintf([cID sprintf('Current sample size: %d, iteration: %.3e\n', tSamp, time)]); 
            end % checkpoint
        end
        fprintf([cID sprintf('Sampling complete, total iteration: %.3e\n', time)]);
        Mx(:,tt,ff,1) = mx;
        My(:,tt,ff,1) = my;
        Mz(:,tt,ff,1) = mz;
        spin_e(:,:,tt,ff) = eSpinT;
        spin_n(:,:,tt,ff) = nSpinT;
        Enuc(:,tt,ff) = Enn;
    end
end
end