function [coef, meas] = MC_sample_GPU(const, ion, params, E_si, E_int, hamI, basis, wav, eSpinT, nSpinT)
% Initialize magnetic moment accumulators on the GPU
Mx = gpuArray.zeros(params.sampSize, length(params.temp), size(params.field, 2));
My = gpuArray.zeros(params.sampSize, length(params.temp), size(params.field, 2));
Mz = gpuArray.zeros(params.sampSize, length(params.temp), size(params.field, 2));
for Sz = 1:params.sampSize
    % Sampling loop
    for tt = 1:length(params.temp)
        beta = gpuArray(1 / (const.J2meV * const.kB * params.temp(tt)));
        for ff = 1:size(params.field, 2)
            field = gpuArray(params.field(:,ff));
            % Perform sampling over the defined interval
            for tInt = 1:params.mIntv
                [E_si(:,tt,ff), E_int(:,tt,ff), ~, wav(:,:,tt,ff), eSpinT(:,:,tt,ff), ~] = thermalize_GPU(const, ion, params,...
                    E_si(:,tt,ff), E_int(:,tt,ff), hamI(:,:,tt,ff), wav(:,:,tt,ff), basis(:,:,tt,ff), eSpinT(:,:,tt,ff), nSpinT(:,:,tt,ff), beta);
                % Update nuclear spins conditionally
                if params.hyp
                    nSpinT(:, :, tt, ff) = therm_nuc_GPU(const, ion, params, beta, field, eSpinT(:, :, tt, ff), nSpinT(:, :, tt, ff));
                end
            end
            % Calculate and accumulate magnetic moments
            Mx(Sz, tt, ff) = rms(eSpinT(:, 1, tt, ff));
            My(Sz, tt, ff) = rms(eSpinT(:, 2, tt, ff));
            Mz(Sz, tt, ff) = rms(eSpinT(:, 3, tt, ff));
        end
    end
    if params.saveOpt == true && mod(Sz, params.sIntv) == 0
        % gather results to CPU to save to file
        coef = gather(wav);
        [meas.Mx, meas.My, meas.Mz] = gather(Mx, My, Mz);
        meas.eSpin = gather(eSpinT);
        meas.nSpin = gather(nSpinT);
        meas.Esi = gather(E_si);
        meas.Eint = gather(E_int);
        meas.Em = sum(E_si + E_int, 1) / size(params.pos,1); % <E> (meV)
        meas.E2m = sum( (E_si + E_int).^2, 1) / size(params.pos,1); % <E^2> (meV^2)
        save(params.DataFile, 'coef', 'meas', '-append');
        fprintf('Data saved!\n');

        % return data back to GPU
        eSpinT = gpuArray(eSpinT);
        nSpinT = gpuArray(nSpinT);
        E_si = gpuArray(E_si);
        E_int = gpuArray(E_int);
        wav = gpuArray(wav);
        Mx = gpuArray(Mx);
        My = gpuArray(My);
        Mz = gpuArray(Mz);
    end
end
coef = gather(wav);
[meas.Mx, meas.My, meas.Mz] = gather(Mx, My, Mz); % gather results to CPU
meas.eSpin = gather(eSpinT);
meas.nSpin = gather(nSpinT);
meas.Esi = gather(E_si);
meas.Eint = gather(E_int);
meas.Em = sum(E_si + E_int, 1) / size(params.pos,1); % <E> (meV)
meas.E2m = sum( (E_si + E_int).^2, 1) / size(params.pos,1); % <E^2> (meV^2)
end