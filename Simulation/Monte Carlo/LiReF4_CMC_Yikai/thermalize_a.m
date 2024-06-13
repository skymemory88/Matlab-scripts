function [dE, accpRate, Esi, Eint, coef, eSpin] = thermalize_a(ion, const, params, beta, Esi, Eint, hamI, coef, basis, eSpin, nSpin)
pos = params.pos;
% site = 1:size(pos,1); % uniform update of the whole lattice
site = randperm(size(pos,1)); % random update of the whole lattice

% generate a random rotation on the Bloch sphere
coords_new = randSph(length(site),'Cartesian');
[alp_new, bet_new, ~] = cart2sph(coords_new(:,1),coords_new(:,2),coords_new(:,3));
alp_new = alp_new + pi; % [-pi pi] -> [0 2pi]
bet_new = bet_new + pi/2; % [-pi/2 pi/2] -> [0 pi]
coef_new = [cos(bet_new'/2) ; sin(bet_new'/2).*exp(1i*alp_new')];

counter = 1; % iteration tracker
change = zeros(length(site),1); % Mark acceptances
dEi = zeros(length(site),1); % local energy change
for ii = site
    % compute new effective spin moment
    wav = squeeze(coef_new(:, ii));
    spx = real(wav' * basis' * ion.Jx * basis * wav);
    spy = real(wav' * basis' * ion.Jy * basis * wav);
    spz = real(wav' * basis' * ion.Jz * basis * wav);
    spin_new = [spx spy spz];
    
    % compute change of single-ion energy
    Esi_new = real(wav' * hamI * wav); % new single-ion energy 
    if params.hyp && ismember(ii, params.isoIdx) % check for hyperfine interaction option
        Esi_new = Esi_new + ion.A(ion.idx) * dot(nSpin(ii,:), spin_new); % include hyperfine energy (meV)
    end
    dE_si = Esi_new - Esi(ii); % single-ion energy

    Eint_new = CntrDip(const.gfac, pos, eSpin, ii, spin_new); % dipole interaction energy change
    if ion.ex(ion.idx) % check exchange interaction strength        
        Eint_new = Eint_new + exchange(params, ion, eSpin, ii, spin_new); % include exchange interaction energy
    end
    dE_int = Eint_new - Eint(ii); % interaction energy change
    
    % total energy difference
    dEi(counter) = dE_si + dE_int;

    % Glauber algorithm
    crit = rand; % update criteria
    prob = 1 / ( 1 + exp( dEi(counter) * beta )); % probability critirion
    if prob >= crit
        Esi(ii) = Esi_new;
        coef(:,ii) = wav;
        eSpin(ii, :) = spin_new;
        Eint(ii) = Eint_new;
        change(counter) = 1;
    end
    counter = counter + 1;
end
% stop thermalization pending on the acceptance rate or the global energy change
dE = sum(dEi .* change); % global energy change
accpRate = sum(change)/length(change); % acceptance rate
end