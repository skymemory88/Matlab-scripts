function [En, wav] = dipole_test(mion, q, range)

renorm = [1, 1, 1]; % renormalization factor
alpha = 0;
Option.demag = false; % demagnetization option
Option.symmetry = true; % use symmetry to infer spin configuration inside the unit cell

theta = 0;
phi = 0;

%      Er       Ho      Yb      Tm      Gd      Y
ion.J = [15/2;    8;      7/2;    6;      7/2;    1];
ion.L = [6;       6;      3;      5;      0;      1];
ion.S = [3/2;     2;      1/2;    1;      7/2;    1];

%Ions' lattice parameters
lattc = [{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Er
           {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      % Ho
%            {[5.175 0 0; 0 5.175 0; 0 0 5.175]}      % Ho test
           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      % Yb
           {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      % Tm
           {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Gd
%            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}     % Gd
           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    % Y

switch mion
    case 'Er'
        eIdx = 1; % element index (1-5): Er Ho Yb Tm Gd Y
        mom = [1 0 0; -1 0 0; -1 0 0; 1 0 0]; % spin moment config
    case 'Ho'
        eIdx = 2;
%         mom = [0 0 1;  0 0 1;  0 0 1; 0 0 1];
        Jxi = cos(phi)*cos(theta);
        Jyi = sin(phi)*cos(theta);
        Jzi = sin(theta);
        mom = [Jxi  Jyi  Jzi
               Jxi  Jyi  Jzi
               Jxi  Jyi  Jzi
               Jxi  Jyi  Jzi];
    case 'Yb'
        eIdx = 3;
        mom = [1 0 0; -1 0 0; -1 0 0; 1 0 0];
    case 'Tm'
        eIdx = 4;
        mom = [1 0 0; -1 0 0; -1 0 0; 1 0 0];
    case 'Gd'
        eIdx = 5;
        mom = [1 0 0; 1 0 0; -1 0 0; -1 0 0];
    case 'Y'
        eIdx = 6;
        mom = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
end
muB = 9.274e-24; % Bohr magneton [J/T]
mu0 = 4e-7*pi; % Vacuum permeability [H/m]
J2meV = 6.24151e+21; % [mev/J]
gfac = mu0*muB^2/4/pi*J2meV*10^30;

lattc = lattc{eIdx};
gL = gLande(ion.L(eIdx),ion.S(eIdx));
Vc = sum(lattc(1,:).*cross(lattc(2,:),lattc(3,:)))*10^-30; % Volume of unit cell [A^-3]

% Initiate J operators
J = ion.J(eIdx);
Jz = diag(J:-1:-J);
Jp = diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm = Jp';
Jx = (Jp+Jm)/2;
Jy = (Jp-Jm)/2i;

% demagnetization parameters
eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
dmag_fac = 4/Vc * mu0/4/pi * muB^2 * J2meV; % demag_fac = N/V * mu_0/4pi * (mu_b)^2
demagn = zeros(3,3,4,4);
demagn_t = ellipsoid_demagn(alpha);
% demagn_t = [7.55 0 0; 0 1.68 0; 0 0 -0.32]; % Thesis "Ultra-low temperature dilatometry" by John L. Dunn
demagn(1,1,:,:) = demagn_t(1,1);
demagn(2,2,:,:) = demagn_t(2,2);
demagn(3,3,:,:) = demagn_t(3,3);

En = double.empty(2*J+1, length(range), 0);
wav = double.empty(2*J+1, 2*J+1, length(range), 0);
parfor ii = 1:length(range)
    worker = getCurrentTask();
    fprintf('Calculation for range: %1$u cells. Core %2$u.\n', range(ii), worker.ID);
    d_sum = dipole_direct(q, range(ii), lattc); % dipole sum
    dipE = gL*(gfac*d_sum + dmag_fac*pi*(eins/3 - Option.demag*demagn));
    for ionn=1:4 % iteracte through each of the four Ho3+ ion in the unit cell
        h_dipol = zeros(1,3);
        for ionm=1:4
            h_dipol = h_dipol + mom(ionm,:)*diag(renorm)*dipE(:,:,ionm,ionn)'; % dipole energy [meV]
        end
        h_dip = -gL*(h_dipol(1)*Jx + h_dipol(2)*Jy + h_dipol(3)*Jz);
        [wv, ee] = eig(h_dip); % [meV]
        ee = real(diag(ee));
        ee = ee-min(ee); % normalize the eigenenergy to the ground state
        [En(:,ii,1), n] = sort(ee);
        wav(:,:,ii,1) = wv(:,n); % reorder the eigenfunctions according to eigenstates
        if Option.symmetry == true; break; end
    end
end

end