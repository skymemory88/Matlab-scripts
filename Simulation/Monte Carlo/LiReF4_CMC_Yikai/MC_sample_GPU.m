function [coef, meas] = MC_sample_GPU(const, ion, params, E_si, hamI, basis, wav, init)
% Initialize magnetic moment accumulators on the GPU
eSpin = init.eSpinT; % intermediate electronic spin configuration after thermalization
nSpin = init.nSpinT; % intermediate nuclear spin configuration after thermalization
pos = params.pos;
pos = gpuArray(pos);

% Sampling loop
for Sz = 1:params.sampSize
    % transfer data to GPU
    E_si = gpuArray(E_si);
    coef = gpuArray(wav);
    eSpinT = gpuArray(eSpin);
    nSpinT = gpuArray(nSpin);

    % Container for magnetic moments along x, y, z
    Mx = gpuArray.zeros(length(params.temp), size(params.field, 2));
    My = gpuArray.zeros(length(params.temp), size(params.field, 2));
    Mz = gpuArray.zeros(length(params.temp), size(params.field, 2));

    for tt = 1:length(params.temp)
        beta = gpuArray(1 / (const.J2meV * const.kB * params.temp(tt)));
        for ff = 1:size(params.field, 2)
            field = gpuArray(params.field(:,ff));
            % Perform sampling over the defined interval
            for tInt = 1:params.mIntv
                [E_si(:,tt,ff), ~, coef(:,:,tt,ff), eSpinT(:,:,tt,ff), ~] = thermalize_GPU(const, ion, params,...
                    E_si(:,tt,ff), hamI(:,:,tt,ff), coef(:,:,tt,ff), basis(:,:,tt,ff), eSpinT(:,:,tt,ff), nSpinT(:,:,tt,ff), beta);
                % Update nuclear spins conditionally
                if params.hyp
                    nSpinT(:, :, tt, ff) = therm_nuc_GPU(const, ion, params, beta, field, eSpinT(:, :, tt, ff), nSpinT(:, :, tt, ff));
                end
            end
            % Calculate and accumulate magnetic moments
            Mx(tt, ff) = rms(eSpinT(:, 1, tt, ff));
            My(tt, ff) = rms(eSpinT(:, 2, tt, ff));
            Mz(tt, ff) = rms(eSpinT(:, 3, tt, ff));
            E_int = dipSum_GPU(const.gfac, pos, eSpinT(:,:,tt,ff));
        end
    end
    % gather results to CPU to save to file
    [wav, E_si] = gather(coef, E_si);
    [Mx, My, Mz] = gather(Mx, My, Mz);
    eSpin = gather(eSpinT);
    nSpin = gather(nSpinT);
    Et_si = reshape(sum(E_si,1), length(params.temp), size(params.field,2));
    E_avg = (Et_si + E_int) / size(params.pos,1);
    if Sz == 1
        meas.Mx = Mx;
        meas.My = My;
        meas.Mz = Mz;
        meas.Esi = E_si;
        meas.Em = E_avg; % <E> (meV)
    elseif Sz > 1
        meas.Mx = cat(3, meas.Mx, Mx);
        meas.My = cat(3, meas.My, My);
        meas.Mz = cat(3, meas.Mz, Mz);
        meas.Esi = cat(3, meas.Esi, E_si);
        meas.Em = cat(3, meas.Em, E_avg); % <E> (meV)
    end
    meas.eSpin = eSpin;
    meas.nSpin = nSpin;
    if params.saveOpt == true && ismember(Sz, 1 + (0:params.sampSize-1)*params.sIntv)
        save(params.DataFile, 'wav', 'meas', '-append');
        fprintf('Data saved!\n');
    end
end
end