function [dE, accpRate, Esi, Edip, coef, spin] = thermalize_a(ion, const, params, tt, Esi, Edip, hamI, coef, basis, spin, spin_n)
pos = params.pos;
% site = 1:size(pos,1); % uniform update of the whole lattice
site = randperm(size(pos,1)); % random update of the whole lattice

% rotations on the Bloch sphere
coords_new = randSph(length(site),'Cartesian');
[alp_new, bet_new, ~] = cart2sph(coords_new(:,1),coords_new(:,2),coords_new(:,3));
alp_new = alp_new + pi; % [-pi pi] -> [0 2pi]
bet_new = bet_new + pi/2; % [-pi/2 pi/2] -> [0 pi]
coef_new = [cos(bet_new'/2) ; sin(bet_new'/2).*exp(1i*alp_new')];

change = zeros(length(site),1); % Mark acceptances
dEi = zeros(length(site),1); % local energy change
for ii = 1:length(site)
    % compute new effective spin moment
    wav = squeeze(coef_new(:,ii));
    spx = real(wav' * basis' * ion.Jx * basis * wav);
    spy = real(wav' * basis' * ion.Jy * basis * wav);
    spz = real(wav' * basis' * ion.Jz * basis * wav);
    spin_new = [spx spy spz];
    
    % compute change of single-ion energy
    Esi_new = real(wav' * hamI * wav); % new single-ion energy 
    % thermalize nuclear spins to the original configuration
    if params.hyp
        Esi_new = Esi_new + ion.A(ion.idx) * sum(spin_n(site(ii),:) .* spin_new); % include hyperfine energy (meV)
    end
    dEs = Esi_new - Esi(site(ii)); % single-ion energy

    % compute dipole interaction energy change
    Edip_new = CntrDip(const, ion, pos, spin, site(ii), spin_new);
    dEd = Edip_new - Edip(site(ii)); % dipole interaction energy
    
    % total energy difference
    dEi(ii) = dEs + dEd;

    % Glauber algorithm
    prob = 0;
    crit = rand; % update criteria
    beta = 1 / (const.J2meV * const.kB * params.temp(tt)); % [meV]
    if params.temp(tt) > 0
        prob = 1 / ( 1 + exp( dEi(ii) * beta )); % probability critirion
    end

    if prob >= crit
        Esi(site(ii)) = Esi_new;
        Edip(site(ii)) = Edip_new;
        coef(:,site(ii)) = wav;
        spin(site(ii), :) = spin_new;
        change(ii) = change(ii) + 1;
    end
end
% stop thermalization pending on the acceptance rate or the global energy change
dE = sum(dEi .* change); % global energy change
accpRate = sum(change)/size(change,1); % acceptance rate
end