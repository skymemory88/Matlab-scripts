function [dE, accpRate, en0, coef, spin] = thermalize(ion, const, params, temp, en0, hamI, coef, basis, spin)
pos = params.pos;
site = 1:size(pos,1); % uniform update of the whole lattice
% site = randperm(size(pos,1)); % random update of the whole lattice

% rotate on the Bloch sphere
coords_new = randSph(length(site),'Cartesian');
[alp_new, bet_new, ~] = cart2sph(coords_new(:,1),coords_new(:,2),coords_new(:,3));
alp_new = alp_new + pi; % [-pi pi] -> [0 2pi]
bet_new = bet_new + pi/2; % [-pi/2 pi/2] -> [0 pi]
coef_new = [cos(bet_new'/2) ; sin(bet_new'/2).*exp(1i*alp_new')];

change = zeros(length(site),1); % Mark acceptances
dEi = zeros(length(site),1); % local energy change
dipE0 = double.empty(length(site),0); % current dipole-dipole interaction energy
for ii = 1:length(site)
    % compute new effective spin moment
    wav = squeeze(coef_new(:,ii));
    spx = real(wav' * basis' * ion.Jx * basis * wav);
    spy = real(wav' * basis' * ion.Jy * basis * wav);
    spz = real(wav' * basis' * ion.Jz * basis * wav);
    spin_new = [spx spy spz];
    
    % compute change of single-ion energy
    En = real(wav' * hamI * wav); % new single-ion energy 
    dEs = En - en0(ii); % single-ion energy

    % compute dipole interaction energy change
    dipE0(ii,1) = CntrDip(const, ion, pos, spin, site(ii), spin(site(ii),:)); % initial dipolar interaction energy
    dEd = CntrDip(const, ion, pos, spin, site(ii), spin_new) - dipE0(ii); % dipole interaction energy
    
    % total energy difference
    dEi(ii) = dEs + dEd;

    % Glauber algorithm
    prob = 0;
    crit = rand; % update criteria
    if temp > 0
        prob = 1 / ( 1 + exp(2*dEi(ii)/const.kB/const.J2meV/temp) ); % probability critirion
    end

    if dEi(ii) < 0 || prob >= crit
        en0(ii) = En;
        coef(:,ii) = wav;
        spin(site(ii), :) = spin_new;
        change(ii) = change(ii) + 1;
    end
end
% stop thermalization pending on the acceptance rate or the global energy change
dE = sum(dEi .* change); % global energy change
accpRate = sum(change)/size(change,1); % acceptance rate
end