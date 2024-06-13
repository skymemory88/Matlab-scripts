function [E_si, E_int, E_tot, eSpin0, nSpin0, params] = initialize_GPU(const, ion, params, coef, basis, hamI)
nnD = min(nonzeros(vecnorm(ion.tau * ion.abc{ion.idx}, 2, 2)));  % Nearest neighbour distance for exchange interaction
E_si = gpuArray(double.empty(size(params.pos,1), length(params.temp), size(params.field, 2), 0)); % single-ion energy
E_int = gpuArray(double.empty(size(params.pos,1), length(params.temp), size(params.field, 2), 0)); % single-ion energy
eSpin0 = zeros(size(params.pos,1), 3, length(params.temp), size(params.field, 2), 'gpuArray'); % initial electornic spin config
nSpin0 = zeros(size(params.pos,1), 3, length(params.temp), size(params.field, 2), 'gpuArray'); % nuclear spin config
E_tot = zeros(params.tEQ, length(params.temp), size(params.field, 2), 'gpuArray'); % total energy
for tt = 1:length(params.temp)
    for ff = 1:size(params.field,2)
        wav = squeeze(coef(:,:,tt,ff));
        bas = squeeze(basis(:,:,tt,ff)); % single-ion wavefunction wav, basis, ion, hamI
        eSpin_temp = arrayfun(@(ii) newSpin_GPU(ion, bas, wav(:,ii), 'electron'), 1:size(params.pos,1), 'UniformOutput', false);
        eSpin0(:,:,tt,ff) = cat(1, eSpin_temp{:}); % concatenate and flatten the results
        E_si(:,tt,ff,1) = real(arrayfun(@(ii) wav(:,ii)' * squeeze(hamI(:,:,tt,ff)) * wav(:,ii), 1:size(params.pos,1))); % initial single-ion energy
        E_int(:,tt,ff, 1) = arrayfun(@(ii) CntrDip_GPU(const.gfac, params.pos, eSpin0(:,:,tt,ff), ii, eSpin0(ii,:,tt,ff)), 1:size(params.pos,1)); % Initial dipole interaction energy
        if ion.ex(ion.idx) % check exchange interaction strength
            E_ex = arrayfun(@(ii) exchange_GPU(params.pos, ion, nnD, eSpin0(:,:,tt,ff), ii, eSpin0(ii,:,tt,ff)), 1:size(params.pos,1), 'UniformOutput', false); % initial exchange interaction energy
            E_ex = cat(1, E_ex{:});
            E_int(:,tt,ff) = E_int(:,tt,ff) + E_ex;
        end
        % Initialize nuclear spins (based on the isotope abundance)
        if params.hyp
            % Number of positions with nuclear moments
            numIso = ceil(size(params.pos, 1) * ion.hyp(ion.idx));
            % Generate random indices for sites with nuclear moments
            isoIdx = randperm(size(params.pos, 1), numIso);
            params.isoIdx = gpuArray(isoIdx);
            eSpin_iso = eSpin0(params.isoIdx, :, tt, ff); % select the isotope sites

            % Nuclear Zeeman interaction matrix for all selected positions
            HzI = -ion.nLande(ion.idx) * const.muN * const.J2meV * (params.field(1, ff)...
                * ion.Ix + params.field(2, ff) * ion.Iy + params.field(3, ff) * ion.Iz);
            ham_temp = arrayfun(@(ii) ion.A(ion.idx) * (eSpin_iso(ii,1) * ion.Ix + eSpin_iso(ii,2)...
                * ion.Iy + eSpin_iso(ii,3) * ion.Iz), 1:numIso, 'UniformOutput', false);
            H_hyp = cat(3, ham_temp{:}); % concatenate all the results
            ham_nuc = H_hyp + repmat(HzI, 1, 1, numIso);  % Expand HzI to match H_hyp dimensions

            % Diagonalization and spin configuration setup
            % NOTE: This step might still need looping unless a batched eig solver that works on GPU exists
            for ii = 1:length(params.isoIdx)
                [eigen_n, ~] = eig(ham_nuc(:, :, ii));
                wav_n = randSphN(1, size(ham_nuc, 1)); % random nuclear spin configuration
                nSpx = real(wav_n' * eigen_n' * ion.Ix * eigen_n * wav_n);
                nSpy = real(wav_n' * eigen_n' * ion.Iy * eigen_n * wav_n);
                nSpz = real(wav_n' * eigen_n' * ion.Iz * eigen_n * wav_n);
                nSpin0(params.isoIdx(ii), :, tt, ff) = [nSpx nSpy nSpz];
                E_si(ii, tt, ff, 1) = E_si(ii, tt, ff, 1) + ion.A(ion.idx)...
                    * dot(nSpin0(params.isoIdx(ii), :, tt, ff), eSpin_iso(ii, :));
            end
        end
    end
end
fprintf('Initialization compelte!\n')
end