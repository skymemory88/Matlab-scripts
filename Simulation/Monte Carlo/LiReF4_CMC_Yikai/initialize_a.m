function [E_si, eSpin0, nSpin0, params] = initialize_a(const, ion, params, coef, basis, hamI)
E_si = double.empty(size(params.pos,1), length(params.temp), size(params.field, 2), 0); % single-ion energy
eSpin0 = zeros(size(params.pos,1), 3, length(params.temp), size(params.field, 2)); % initial electornic spin config
nSpin0 = zeros(size(params.pos,1), 3, length(params.temp), size(params.field, 2)); % nuclear spin config
for tt = 1:length(params.temp)
    for ff = 1:size(params.field,2)
        bas = squeeze(basis(:,:,tt,ff)); % single-ion wavefunction
        for ii = 1:size(params.pos,1)
            % initialize electronic spins
            wav = squeeze(coef(:,ii,tt,ff)); % Unitary rotation
            spx = real(wav' * bas' * ion.Jx * bas * wav)';
            spy = real(wav' * bas' * ion.Jy * bas * wav)';
            spz = real(wav' * bas' * ion.Jz * bas * wav)';
            eSpin0(ii,:,tt,ff) = [spx, spy, spz]; % initial electronic spin config
            E_si(ii,tt,ff,1) = real(wav' * squeeze(hamI(:,:,tt,ff)) * wav); % initial single-ion energy
        end
        % Initialize nuclear spins (based on the isotope abundance)
        if params.hyp
            numIso = ceil(size(params.pos, 1) * ion.hyp(ion.idx));
            params.isoIdx = randperm(size(params.pos,1), numIso); % sites with nuclear moments
            HzI = const.Ngfac * (params.field(1, ff)*ion.Ix + params.field(2, ff)*ion.Iy + params.field(3, ff)*ion.Iz); % nuclear Zeeman interaction (meV)
            for jj = params.isoIdx
                spx = eSpin0(jj,1,tt,ff);
                spy = eSpin0(jj,2,tt,ff);
                spz = eSpin0(jj,3,tt,ff);
                H_hyp = ion.A(ion.idx) * (spx*ion.Ix + spy*ion.Iy + spz*ion.Iz); % hyperfine interaction (meV)
                ham_nuc = H_hyp + HzI; % total nuclear spin hamiltonian
                [eigen_n, ~] = eig(ham_nuc); % diagonalize the nuclear hamiltonian
                wav_n = randSphN(1, size(ham_nuc,1)); % random point of the nuclear spin configuration
                spx_n = real(wav_n' * eigen_n' * ion.Ix * eigen_n * wav_n)';
                spy_n = real(wav_n' * eigen_n' * ion.Iy * eigen_n * wav_n)';
                spz_n = real(wav_n' * eigen_n' * ion.Iz * eigen_n * wav_n)';
                nSpin0(jj,:,tt,ff) = [spx_n spy_n spz_n]; % initial nuclear spin config
                E_si(jj,tt,ff,1) = E_si(jj,tt,ff,1) + ion.A(ion.idx) *...
                    dot(nSpin0(jj,:,tt,ff), eSpin0(jj,:,tt,ff,1)); % include hyperfine interaction
            end
        end
    end
end
fprintf('Initialization compelte!\n')
end