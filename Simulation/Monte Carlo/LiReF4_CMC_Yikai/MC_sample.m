function [Mx, My, Mz, coef, spin, en0] = MC_sample(const, ion, params, en0, hamI, basis, wv, spin0)
pos = params.pos;

spin = double.empty(size(pos,1), 3, length(params.temp), size(params.field, 2), 0); % final spin config
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
        wav = squeeze(wv(:,:,tt,ff)); % single-ion hamiltonian eigenvectors
        spinT = squeeze(spin0(:,:,tt,ff)); % spin configuration
        bas = squeeze(basis(:,:,tt,ff)); % Ising basis

        time = 1; % calculation time tracker
        tSamp = 1; % sample size tracker
        mx = double.empty(params.sIntv,0);
        my = double.empty(params.sIntv,0);
        mz = double.empty(params.sIntv,0);
        while tSamp < params.sIntv
            for tInt = 1:params.mIntv
                [~, ~, en0(:,tt,ff), wav, spinT] = thermalize_a(ion, const, params, params.temp(tt), en0(:,tt,ff), hamI(:,:,tt,ff), wav, bas, spinT);
                time = time + 1;
            end
            % Measurement
            mx(tSamp,1) = mean(spinT(:,1),1); % <Jx>
            my(tSamp,1) = mean(spinT(:,2),1); % <Jy>
            mz(tSamp,1) = mean(spinT(:,3),1); % <Jz>
            tSamp = tSamp + 1;

            if mod(time, 20) == 0; fprintf([cID sprintf('Current sample size: %d, iteration: %.3e\n', tSamp, time)]); end % checkpoint
        end
        fprintf([cID sprintf('Sampling complete, total iteration: %.3e\n', time)]);
        Mx(:,tt,ff,1) = mx;
        My(:,tt,ff,1) = my;
        Mz(:,tt,ff,1) = mz;
        coef(:,:,tt,ff) = wav;
        spin(:,:,tt,ff,1) = spinT;
    end
end
end