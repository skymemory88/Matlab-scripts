function [ion, eee, dE, Js, Is, Jqs, fff, ttt, vvv] = gen_Ising(ttt, fff, theta, phi, plotopt, saveopt)
clearvars -except ttt fff theta phi plotopt saveopt
format long g;

Options.plot = 1; % 1. 3D scatter plot; 2. 2D color map.
Options.pProm = 1e-2; % minimum peak prominance for phase boundary search
Options.intType = 'exchange'; % interaction type: 1. 'exchange' (short range). 2. 'dipole' (long range)
    ion.ex = [-1]; % exchange interaction [meV]
    ion.dpRng = 100; % dipole summation range
Options.damp = 0.0; % damping factor for the iterative update
Options.hyperF = true; % Hyperfine interaction option
    ion.hyp = -[0.1 0.1 1.0]; % Hyperfine interaction [meV]
Options.RPA = true; % random phase approximation
    qz = linspace(0, 1, 200); % q-vector for RPA calculations (0: BZ center, 1: BZ boundary)
    qvec = [zeros(size(qz))' zeros(size(qz))' qz'];
Options.plot = plotopt;
Options.save = saveopt;
Options.filepath = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing',...
        '\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities',...
        '\toy_model'];
if size(fff,1) > length(ttt)
    Options.scanMode = 'field';
else
    Options.scanMode = 'temp';
end

theta = theta*pi/180; % convert to radian, deviation angle from normal vector
phi = phi*pi/180; % convert to radian, in-plane angle from a/b axis

Bx = fff*cos(phi)*cos(theta);
By = fff*sin(phi)*cos(theta);
Bz = fff*sin(theta);
fff = [Bx' By' Bz'];

ion.prop = 1; % isotope proportion
ion.J = 1/2;
ion.gLande = 1; % Lande factor
ion.I = 1/2;
ion.nLande = 0; % nuclear Lande factor

sigx = [0   1
        1   0];
sigy = [ 0   -1i
        1i  0];
sigz = [1   0
        0  -1];

Sx = ion.J * sigx;
Sy = ion.J * sigy;
Sz = ion.J * sigz;

if Options.hyperF == true   
    Ix = ion.I * sigx;
    Iy = ion.I * sigy;
    Iz = ion.I * sigz;

    Sx = kron(Sx, eye(size(Ix,1)));
    Sy = kron(Sy, eye(size(Iy,1)));
    Sz = kron(Sz, eye(size(Iz,1)));

    Ix = kron(eye(size(sigx,1)), Ix);
    Iy = kron(eye(size(sigy,1)), Iy);
    Iz = kron(eye(size(sigz,1)), Iz);
else
    Ix = 0;
    Iy = 0;
    Iz = 0;
    ion.I = 0;
end

% lattice constants [Ang]
ion.abc = [10  0   0;
           0  10   0;
           0   0  10];
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

% Interaction dimension regulation (isotropic, diagonal, Ising)
% ion.mat = [1  1  1;
%            1  1  1;
%            1  1  1]; % isotropic
% ion.mat = [1  0  0;
%            0  1  0;
%            0  0  1]; % diagonal
ion.mat = [0  0  0;
           0  0  0;
           0  0  1]; % Ising
       
kB = 8.61733e-2; % Boltzmann constant [meV/K]
Jq = zeros(3,3,size(qvec,1));
switch Options.intType
    case 'exchange'
        eMeas = abs(max(ion.ex)); % reference energy
%         Jij = diag([0  0  ion.ex]); % Ising chain
        [Jij, ~] = MF_exchange([0 0 0], ion.ex, ion.abc, ion.tau); % exchange coupling sum
        Jij = Jij .* ion.mat; % truncate interaction matrix
        ion.int = 'exchange';
        if Options.RPA == true
            for nq = 1:size(qvec,1)
                [Jqt, ~] = MF_exchange(qvec(nq,:), ion.ex, ion.abc, ion.tau); % exchange coupling sum
                Jqt = Jqt .* ion.mat; % truncate interaction matrix
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
        Jij = gfac * 1e30 * Jij .* ion.mat; % truncate interaction matrix
        eMeas = abs(Jij(3,3)); % reference energy
        ion.int = 'dipole';
        if Options.RPA == true
            for nq = 1:size(qvec,1)
                [Jqt, ~] = MF_dipole(qvec(nq,:), ion.dpRng, ion.abc, ion.tau); % exchange coupling sum
                Jqt = gfac * 1e30 * Jqt .* ion.mat; % truncate interaction matrix
                for ionn = 1:size(ion.tau,1)
                    for ionm = 1:size(ion.tau,1)
                        Jq(:,:,nq) = Jq(:,:,nq) + Jqt(:,:,ionn,ionm);
                    end
                end
            end
        end

end

if Options.hyperF == true
    A = eMeas * ion.hyp;
else
    A = 0;
end

delta = 1e-7; % convergence criteria
iterator = 1; % linear iterator
Js = zeros(size(fff,1), length(ttt), 3);
Is = zeros(size(fff,1), length(ttt), 3);
Jqs = zeros(size(fff,1), length(ttt), 1);
vector = zeros(size(fff,1)*length(ttt), 3);
map = zeros(size(fff,1), length(ttt));
eee = double.empty(size(fff,1), length(ttt), size(Sz,2), 0);
vvv = double.empty(size(fff,1), length(ttt), size(Sz,2), size(Sz,2), 0);
for ii = 1:size(fff,1)
    for jj = 1:length(ttt)
        S_mf = [0  0  ion.J]; % initial guesses for spin operators
        beta = 1/(kB*ttt(jj));
        counter = 1; % while loop iterator
        En_old = 0;
        while true
            h_int = ion.gLande^2 * S_mf * Jij; % spin-spin interaction part 1
            if any(A) % Case: finite hyperfine interaction
                Hint = h_int(1)*Sx + h_int(2)*Sy + h_int(3)*Sz; % spin-spin interaction
                Hzee =  -ion.gLande * (Bx(ii)*Sx + By(ii)*Sy + Bz(ii)*Sz); % Zeeman interaction
                HzeeI = -ion.nLande * (Bx(ii)*Ix + By(ii)*Iy + Bz(ii)*Iz); % nuclear Zeeman interaction
                Hyp = A(1)*Sx*Ix + A(2)*Sy*Iy + A(3)*Sz*Iz; % hyperfine interaction
                ham = Hint + Hzee + HzeeI + Hyp;
                [wav, En] = eig(ham);
                En = real(diag(En));
                [En, n] = sort(En); % sort the energy in ascending order
                wav = wav(:,n); % sort the eigen-vector accordingly
                eee(ii,jj,:,1) = squeeze(En);
                vvv(ii,jj,:,:,1) = squeeze(wav);
                En = En - min(En); % normalize the eigenenergy to the ground state
                Z = exp(-En*beta)/sum(exp(-En*beta)); % Calculate the partition function weight
                if ttt(jj) == 0
                    newAvSx = real( wav(:,1)'*Sx*wav(:,1) ); % Calculate the expectation value of Jz
                    newAvSy = real( wav(:,1)'*Sy*wav(:,1) ); % Calculate the expectation value of Jz
                    newAvSz = real( wav(:,1)'*Sz*wav(:,1) ); % Calculate the expectation value of Jz
                
                    newAvIx = real(wav(:,1)'*Ix*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvIy = real(wav(:,1)'*Iy*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvIz = real(wav(:,1)'*Iz*wav(:,1))'; % Calculate the expectation value of Jz
                else
                    newAvSx = real( diag(wav'*Sx*wav) )'*Z; % Calculate the expectation value of Jz
                    newAvSy = real( diag(wav'*Sy*wav) )'*Z; % Calculate the expectation value of Jz
                    newAvSz = real( diag(wav'*Sz*wav) )'*Z; % Calculate the expectation value of Jz
                
                    newAvIx = real(wav(:,1)'*Ix*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvIy = real(wav(:,1)'*Iy*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvIz = real(wav(:,1)'*Iz*wav(:,1))'; % Calculate the expectation value of Jz
                end
                newAvI = [newAvIx  newAvIy  newAvIz];
            else % Case: no hyperfine interaction
%                 H_ani = 5e-3*fff(ii)* Sz; % anisotropy term
                Hzee = -ion.gLande * Bx(ii)*eMeas * Sx - ion.gLande * By(ii)*eMeas * Sy - ion.gLande * Bz(ii)*eMeas * Sz;
                ham = ion.gLande * (h_int(1)*Sx + h_int(2)*Sy + h_int(3)*Sz) + Hzee;
                [wav, En] = eig(ham);
                En = real(diag(En));
                eee(ii,jj,:,1) = squeeze(En);
                vvv(ii,jj,:,:,1) = squeeze(wav);
                En = En-min(En); % normalize the eigenenergy to the ground state
                Z = exp(-En*beta)/sum(exp(-En*beta)); % Calculate the partition function weight
                if ttt(jj) == 0
                    newAvSx = real(wav(:,1)'*Sx*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvSy = real(wav(:,1)'*Sy*wav(:,1))'; % Calculate the expectation value of Jz
                    newAvSz = real(wav(:,1)'*Sz*wav(:,1))'; % Calculate the expectation value of Jz
                else
                    newAvSx = real(diag(wav'*Sx*wav))'*Z; % Calculate the expectation value of Jz
                    newAvSy = real(diag(wav'*Sy*wav))'*Z; % Calculate the expectation value of Jz
                    newAvSz = real(diag(wav'*Sz*wav))'*Z; % Calculate the expectation value of Jz
                end
            end
            newAvS = [newAvSx  newAvSy  newAvSz];
                
            if ~any( abs(En - En_old) >= delta ) && ~any( abs(newAvS(3) - S_mf(3)) >= delta )
                break
            elseif counter > 1e5
                fprintf('Iteration exceeded limits, dE = %.2e \n', nonzeros(En - En_old));
                break
            else
                En_old = En; % store energy of the last round
                S_mf = Options.damp*S_mf + (1-Options.damp)*newAvS; % update the order parameter (internal field)
                counter = counter + 1;
                if mod(counter, 2e4) == 0
                   fprintf('Iternation: %u, Magnetic field: %.2f, Temperature: %.2f\n', counter, Bx(ii), ttt(jj)) 
                end
            end
        end
        Js(ii,jj,:) = S_mf;
        if exist('newAvI','var'); Is(ii,jj,:) = newAvI; end
        vector(iterator,:) = [ttt(jj), Bx(ii), S_mf(3)]; % Phase diagram with Sz as the order parameter
        map(ii,jj) = S_mf(3); % use Sz as the order parameter
        iterator = iterator + 1;
        %if mod(iterater, 1000) == 0;
        %    iterater, temp, hT
        %end
    end
    if any(ismember(round(linspace(1,size(fff,1),20)),ii))
        fprintf('Calculating for |B/J| = %.2f, progress: %.1f%% \n', Bx(ii), ii/size(fff,1)*100);
    end 
end
Js = squeeze(Js);
Is = squeeze(Is);
eee = squeeze(eee); % eigen-states
ee = eee - min(eee,[],2); % normalize to the ground state energy
dE = zeros(nchoosek(size(ee,2),2),size(fff,1)); % excitation modes
for ii = 1:size(fff,2)
    [en, em] = meshgrid(squeeze(ee(ii,:)),squeeze(ee(ii,:)));
    delta = triu(en-em);
    delta = delta(delta>0);
    dE(:,ii) = padarray(delta, size(dE,1)-size(delta,1), 'pre');
%     dE(:,ii) = delta(:);
end
vvv = squeeze(vvv);

if Options.plot == true
    Bc = zeros(length(ttt),1); % critical fields (at fixed temperatre)
    for ii = 1:length(ttt)
        dif = real(diff(map(:,ii))); % crude differentiation
        [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Options.pProm);
        if isempty(idx)
            Bc(ii) = 0;
        else
            Bc(ii) = Bx(idx);
        end
    end
    TT = ttt(Bc~=0)'; % remove zeros and nonsensical values

    Tc = zeros(size(fff,1),1); % critical temperatures (at fixed field)
    for ii = 1:size(fff,1)
        dif = real(diff(map(ii,:))); % crude differentiation
        [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Options.pProm);
        if isempty(idx)
            Tc(ii) = 0;
        else
            Tc(ii) = ttt(idx);
        end
    end
    BB = Bx(Tc~=0)'; % remove zeros and nonsensical values
    phB = [TT nonzeros(Bc); nonzeros(Tc) BB]./eMeas; % Phase boundary

    figure
    box on
    hold on
    switch Options.plot
        case 1
            scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
            view([0 90]);
        case 2
            temps = repmat(ttt, size(fff,1), 1);
            fields = repmat(Bx, 1, length(ttt));
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
        ttt, max(vecnorm(fff,2,2)), any(ion.hyp)));
    fileobj = fullfile(Options.filepath,char(tit));
    save(fileobj,'ttt','fff','eee','vvv','ion','-v7.3')
end
end