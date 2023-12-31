function [ion, eee, dE, Js, Is, fff, ttt, vvv, addon] = Hsbrg_mod(ttt, fff, theta, phi, plotopt, saveopt)
clearvars -except ttt fff theta phi plotopt saveopt

Options.pProm = 1e-2; % minimum peak prominance for phase boundary search
Options.intType = 'exchange'; % interaction type: 1. 'exchange' (short range). 2. 'dipole' (long range)
    ion.ex = [-1]; % exchange interaction [meV]
    ion.dpRng = 100; % dipole summation range
Options.damp = 0.0; % damping factor for the iterative update
Options.IsingProj = true; % Project to Ising space
Options.hyperF = true; % Hyperfine interaction option
    ion.hyp = -[0.8 0.8 0.8]; % Hyperfine interaction [meV]
Options.RPA = true; % random phase approximation
    omega = linspace(0,30,801);
%     qz = linspace(0, 1, 200); % q-vector for RPA calculations (0: BZ center, 1: BZ boundary)
%     qvec = [zeros(size(qz))' zeros(size(qz))' qz'];
    qvec = [0 0 0];
Options.plot = plotopt; % 0, no plot; 1. 3D scatter plot; 2. 2D color map.
Options.save = saveopt;
Options.filepath = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program',...
    '\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities\toy_model'];
if size(fff,1) > length(ttt)
    Options.scanMode = 'field';
else
    Options.scanMode = 'temp';
end

ion.prop = 1; % isotope proportion
ion.J = 1; % electronic spin length
ion.gLande = 1; % Lande factor
Dz = -11; % Ising anisotropy
ion.I = 1/2; % nuclear spin length
ion.nLande = 0; % nuclear Lande factor

% Define operators:
Sz = diag(ion.J:-1:-ion.J);
Sp = diag(sqrt((ion.J-(ion.J-1:-1:-ion.J)).*(ion.J+1+(ion.J-1:-1:-ion.J))),1);
Sm = Sp';
Sx = (Sp + Sm)/2;
Sy = (Sp - Sm)/2i;

if Options.hyperF == true   
    Iz = diag(ion.I:-1:-ion.I);
    Ip = diag(sqrt((ion.I-(ion.I-1:-1:-ion.I)).*(ion.I+1+(ion.I-1:-1:-ion.I))),1);
    Im = Ip';
    Ix = (Ip + Im)/2;
    Iy = (Ip - Im)/2i;

    Sx = kron(Sx, eye(2*ion.I+1));
    Sy = kron(Sy, eye(2*ion.I+1));
    Sz = kron(Sz, eye(2*ion.I+1));

    Ix = kron(eye(2*ion.J+1), Ix);
    Iy = kron(eye(2*ion.J+1), Iy);
    Iz = kron(eye(2*ion.J+1), Iz);
else
    ion.I = 0;
    Ix = 0;
    Iy = 0;
    Iz = 0;
end

% lattice constants [Ang]
ion.abc = [1   0   0;
           0   0   0;
           0   0   0];
% ion.abc = [5.175    0      0;
%              0    5.175    0;
%              0      0    10.75]; % LiHoF4 lattice

% basis vector within the unit cell
ion.tau = [0   0   0];
% ion.tau = [ 0    0    0;
%            1/2  1/2  1/2]; % bcc
% ion.tau = [ 0     0     0;
%             0    1/2   1/4;
%            1/2   1/2   1/2;
%            1/2    0    3/4]; % LiHoF4 basis

kB = 8.61733e-2; % Boltzmann constant [meV/K]
Jq = zeros(3,3,size(qvec,1));
switch Options.intType
    case 'exchange'
        eMeas = abs(max(ion.ex)); % reference energy
        [Jij, ~] = MF_exchange([0 0 0], ion.ex, ion.abc, ion.tau); % exchange coupling sum
        ion.int = 'exchange';
        if Options.RPA == true
            for nq = 1:size(qvec,1)
                [Jqt, ~] = MF_exchange(qvec(nq,:), ion.ex, ion.abc, ion.tau); % exchange coupling sum
                for ionn = 1:size(ion.tau,1)
                    for ionm = 1:size(ion.tau,1)
                        Jq(:,:,nq) = Jq(:,:,nq) + Jqt(:,:,ionn,ionm);
                    end
                end
            end
        end
    case 'dipole'
        muB = 9.274e-24; % [J/T]
        mu0 = 4e-7*pi; % [H/m]
        J2meV = 6.24151e+21; % [mev/J]
        gfac = muB^2 * mu0 * J2meV / (4*pi); % [meV*Ang^-3] (dipole_direct gives [Ang^3])
        [Jij, ~] = MF_dipole([0 0 0], ion.dpRng, ion.abc, ion.tau); % dipole coupling sum
        Jij = gfac * 1e30 * Jij; % truncate interaction matrix
        eMeas = abs(Jij(3,3)); % reference energy
        ion.int = 'dipole';
        if Options.RPA == true
            for nq = 1:size(qvec,1)
                [Jqt, ~] = MF_dipole(qvec(nq,:), ion.dpRng, ion.abc, ion.tau); % exchange coupling sum
                Jqt = gfac * 1e30 * Jqt; % truncate interaction matrix
                for ionn = 1:size(ion.tau,1)
                    for ionm = 1:size(ion.tau,1)
                        Jq(:,:,nq) = Jq(:,:,nq) + Jqt(:,:,ionn,ionm);
                    end
                end
            end
        end

end
theta = theta*pi/180; % convert to radian, deviation angle from normal vector
phi = phi*pi/180; % convert to radian, in-plane angle from a/b axis

Bx = fff*cos(phi)*cos(theta) * eMeas;
By = fff*sin(phi)*cos(theta) * eMeas;
Bz = fff*sin(theta) * eMeas;
fff = [Bx' By' Bz'] ./ eMeas;

if Options.hyperF == true
    A = eMeas * ion.hyp;
else
    A = 0;
end

delta = 1e-7; % convergence criteria
iterator = 1; % linear iterator
Js = zeros(size(fff,1), length(ttt), 3);
Is = zeros(size(fff,1), length(ttt), 3);
vector = zeros(size(fff,1)*length(ttt), 3);
map = zeros(size(fff,1), length(ttt));

eee = double.empty(size(fff,1), length(ttt), size(Sz,2), 0); % eigen-energies
vvv = double.empty(size(fff,1), length(ttt), size(Sz,2), size(Sz,2), 0); % eigen-functions
BxI = double.empty(size(fff,1), length(ttt), 0); % effective transverse field from Ising projection
dE = zeros(nchoosek(size(Sz,2),2),size(fff,1), length(ttt)); % 1st order excitation modes
deno = double.empty(size(fff,1), length(omega), 0);
poles = zeros(size(fff,1), size(Sz,2)-1);
modes = zeros(size(fff,1), size(Sz,2)-1);
for nb = 1:size(fff,1)
    for nt = 1:length(ttt)
        S_mf = [0  0  ion.J]; % initial guesses for spin operators
        beta = 1/(kB*ttt(nt));
        counter = 1; % while loop iterator
        En_old = 0;
        while true
            h_int = S_mf * Jij; % spin-spin interaction part 1
            if any(A) % Case: finite hyperfine interaction
                Hint = h_int(1)*Sx + h_int(2)*Sy + h_int(3)*Sz; % spin-spin interaction
                Hzee =  -ion.gLande * (Bx(nb)*Sx + By(nb)*Sy + Bz(nb)*Sz); % Zeeman interaction
                HzeeI = -ion.nLande * (Bx(nb)*Ix + By(nb)*Iy + Bz(nb)*Iz); % nuclear Zeeman interaction
                Hyp = A(1)*Sx*Ix + A(2)*Sy*Iy + A(3)*Sz*Iz; % hyperfine interaction
                H_ani = Dz*Sz.^2; % anisotropy term
                ham = Hint + Hzee + HzeeI + Hyp + H_ani; % full hamiltonian
%                 ham = Hzee + HzeeI + Hyp + H_ani; % no spin-spin interaction
                [wav, En] = eig(ham);
                En = squeeze(real(diag(En)));
                [En, n] = sort(En); % sort the energy in ascending order
                wav = wav(:,n); % sort the eigen-vector accordingly
                eee(nb,nt,:,1) = squeeze(En);
                vvv(nb,nt,:,:,1) = squeeze(wav);
                En = En - min(En); % normalize the eigenenergy to the ground state
                Z = exp(-En*beta)/sum(exp(-En*beta)); % Calculate the partition function weight
                if ttt(nt) == 0
                    newAvSx = real( wav(:,1)' *Sx* wav(:,1) ); % Calculate the expectation value of Jz
                    newAvSy = real( wav(:,1)' *Sy* wav(:,1) ); % Calculate the expectation value of Jz
                    newAvSz = real( wav(:,1)' *Sz* wav(:,1) ); % Calculate the expectation value of Jz
                
                    newAvIx = real( wav(:,1)' *Ix* wav(:,1) )'; % Calculate the expectation value of Jz
                    newAvIy = real( wav(:,1)' *Iy* wav(:,1) )'; % Calculate the expectation value of Jz
                    newAvIz = real( wav(:,1)' *Iz* wav(:,1) )'; % Calculate the expectation value of Jz
                else
                    newAvSx = real( diag(wav' *Sx* wav) )'*Z; % Calculate the expectation value of Jz
                    newAvSy = real( diag(wav' *Sy* wav) )'*Z; % Calculate the expectation value of Jz
                    newAvSz = real( diag(wav' *Sz* wav) )'*Z; % Calculate the expectation value of Jz
                
                    newAvIx = real( diag(wav' *Ix* wav) )'*Z; % Calculate the expectation value of Jz
                    newAvIy = real( diag(wav' *Iy* wav) )'*Z; % Calculate the expectation value of Jz
                    newAvIz = real( diag(wav' *Iz* wav) )'*Z; % Calculate the expectation value of Jz
                end
                newAvI = [newAvIx  newAvIy  newAvIz];
            else % Case: no hyperfine interaction
                Hint = h_int(1)*Sx + h_int(2)*Sy + h_int(3)*Sz; % spin-spin interaction
                Hzee = -ion.gLande * (Bx(nb)*Sx + By(nb)*Sy + Bz(nb)*Sz); % Zeeman term
                H_ani = Dz*Sz.^2; % anisotropy term
                ham = Hint + Hzee + H_ani; % full hamiltonian
%                 ham = Hzee + H_ani; % no spin-spin interaction
                [wav, En] = eig(ham);
                En = squeeze(real(diag(En)));
                [En, n] = sort(En); % sort the energy from lowest to the highest
                wav = wav(:,n); % sort the eigen-vectors in its basis accordingly
                eee(nb,nt,:,1) = En;
                vvv(nb,nt,:,:,1) = squeeze(wav);
                En = En - min(En); % normalize the eigenenergy to the ground state
                Z = exp(-En*beta)/sum(exp(-En*beta)); % Calculate the partition function weight
                if ttt(nt) == 0
                    newAvSx = real(wav(:,1)' *Sx* wav(:,1))'; % Calculate the expectation value of Jz
                    newAvSy = real(wav(:,1)' *Sy* wav(:,1))'; % Calculate the expectation value of Jz
                    newAvSz = real(wav(:,1)' *Sz* wav(:,1))'; % Calculate the expectation value of Jz
                else
                    newAvSx = real(diag(wav' *Sx* wav))'*Z; % Calculate the expectation value of Jz
                    newAvSy = real(diag(wav' *Sy* wav))'*Z; % Calculate the expectation value of Jz
                    newAvSz = real(diag(wav' *Sz* wav))'*Z; % Calculate the expectation value of Jz
                end
            end
            newAvS = [newAvSx  newAvSy  newAvSz];
            
            if ~any( abs(En - En_old) >= delta ) && ~any( abs(newAvS(3) - S_mf(3)) >= delta )
                [en, em] = meshgrid(squeeze(En), squeeze(En));
                de = triu(en-em);
                de = unique(de(de>0));
                dE(:,nb,nt) = padarray(de, size(dE,1)-size(de,1), 'pre');
                break
            elseif counter > 1e5              
                [en, em] = meshgrid(squeeze(En), squeeze(En));
                de = triu(en-em);
                de = unique(de(de>0));
                dE(:,nb,nt) = padarray(de, size(dE,1)-size(de,1), 'pre');
                fprintf('Iteration exceeded limits, dE = %.2e \n', nonzeros(En - En_old));
                break
            else
                En_old = En; % store energy of the last round
                S_mf = Options.damp*S_mf + (1-Options.damp)*newAvS; % update the order parameter (internal field)
                counter = counter + 1;
                if mod(counter, 2e4) == 0
                   fprintf('Iternation: %u, Magnetic field: %.2f, Temperature: %.2f\n', counter, Bx(nb)/eMeas, ttt(nt)) 
                end
            end
        end
        if Options.IsingProj == true
            wvI = wav(:,1:2);
            SzI = wvI' * Sz * wvI;
            [uni, ~] = eig(SzI);
            hamI = uni'*wvI' * ham * wvI * uni;
            BxI(nb,nt,1) = 2*mean([abs(hamI(2,1)) abs(hamI(1,2))]);
        end
        Js(nb,nt,:) = S_mf;
        if exist('newAvI','var'); Is(nb,nt,:) = newAvI; end
        vector(iterator,:) = [kB*ttt(nt)/eMeas, Bx(nb), S_mf(3)]; % Phase diagram with Sz as the order parameter
        map(nb,nt) = S_mf(3); % use Sz as the order parameter
        iterator = iterator + 1;
    end
    if any(ismember(round(linspace(1,size(fff,1),20)),nb))
        fprintf('Calculating for |B/J| = %.2f, progress: %.1f%% \n', Bx(nb), nb/size(fff,1)*100);
    end
    if Options.RPA == true
        for nq = 1:size(qvec,1)
%             Smnx = wav' * Sx * wav; % <Sx>_mn
%             Smnx = Smnx(2:end,1); % off diagonal elements between <n| and |0>
%             Smny = wav' * Sy * wav; % <Sy>_mn
%             Smny = Smny(2:end,1);
            Smnz = wav' * Sz * wav; % <Sz>_mn
            Smnz = Smnz(2:end,1);
            Jzz = abs(squeeze(Jq(3,3,nq))); % Ising element of q=0 interaction tensor
%             Jab = abs(squeeze(Jq(:,:,nq))); % Ising element of q=0 interaction tensor
            dE0 = En(2:end); % MF excitation modes at 0 K
            for nw = 1:length(omega) % search for RPA poles
                dE2 = (dE0.^2 - omega(nw)^2);
                deno(nb,nw,1) = prod(dE2)*(1 - (2*dE0') * (abs(Smnz).^2 ./dE2) * Jzz); % Ising interaction
            end
            pks = squeeze(mag2db(1./abs(deno(nb,:,1)))); % search for poles
            [~,loc] = findpeaks(pks, 'SortStr', 'none');
            poles(nb,1:length(loc)) = flip(omega(loc));
            for kk = 2:size(poles,2)
                modes(:,kk-1) = poles(:,kk-1)-poles(:,kk);
            end
        end
    end
end
Js = squeeze(Js);
Is = squeeze(Is);
dE = squeeze(dE);
vvv = squeeze(vvv);
ttt = kB*ttt/eMeas;
if Options.IsingProj == true; addon.Bx = squeeze(BxI);
if Options.RPA ==true
    addon.modes = squeeze(poles);
    addon.deno = squeeze(deno);
    addon.dE = squeeze(modes);
else
    addon = {};
end

if Options.plot ~= 0
    Bc = zeros(length(ttt),1); % critical fields (at fixed temperatre)
    for nb = 1:length(ttt)
        dif = real(diff(map(:,nb))); % crude differentiation
        [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Options.pProm);
        if isempty(idx)
            Bc(nb) = 0;
        else
            Bc(nb) = Bx(idx);
        end
    end
    TT = ttt(Bc~=0)'; % remove zeros and nonsensical values

    Tc = zeros(size(fff,1),1); % critical temperatures (at fixed field)
    for nb = 1:size(fff,1)
        dif = real(diff(map(nb,:))); % crude differentiation
        [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Options.pProm);
        if isempty(idx)
            Tc(nb) = 0;
        else
            Tc(nb) = ttt(idx);
        end
    end
%     phB = [kB*TT./eMeas nonzeros(Bc); kB*nonzeros(Tc)./eMeas Bx(Tc~=0)']; % Phase boundary
    phB = [TT nonzeros(Bc); nonzeros(Tc) Bx(Tc~=0)']; % Phase boundary

    figure
    box on
    hold on
    switch Options.plot
        case 1
            scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
            view([0 90]);
            scatter3(phB(:,1), phB(:,2), ones(size(phB,1),1), '.r');
        case 2
            temps = repmat(kB*ttt/eMeas, size(fff,1), 1);
            fields = repmat(Bx', 1, length(ttt));
            phase = pcolor(temps, fields, map);
            phase.EdgeColor = 'none';
            plot(phB(:,1), phB(:,2), '.r');
    end
    set(gca, 'FontSize', 12)
    xlabel('Temperature (K)')
    ylabel('Magnetic Field (T)')
end

if Options.save == true
    tit = strcat("Hscan_toy_", sprintf('%1$3.3fK_%2$.1fT_hp=%3$.2f.mat',...
        ttt, max(vecnorm(fff,2,2)), Options.hyperF));
    fileobj = fullfile(Options.filepath,char(tit));
    save(fileobj,'ttt','fff','eee','vvv','ion','Js','Is','-v7.3')
end
end