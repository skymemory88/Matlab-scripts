function [init, meas, coef, ion, params] = LiReF4_CMC_Yikai_GPU(mion, dims, iState, hyperfineOpt, loadOpt, plotOpt, saveOpt)
% Argument: mion, loadopt, dim, loadOpt, plotOpt, saveOpt
% mion: Magnetic ion type
% dims: lattice size (1 x 3 array)
% loadOpt: 1. load existing file or, 2. start from scratch
% plotOpt: plotting option: True/False
% saveOpt: saving option: True/False

% set base folder depending on the OS
if ispc
    filepath = ['G:\My Drive\File sharing\PhD research\Li', mion,...
        'F4 project\Data\Simulations\MATLAB\Monte Carlo'];
elseif ismac
    filepath = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
        'PhD program/Research projects/Li', mion, 'F4 project/Data/Simulations/Matlab/Monte Carlo/'];
else
    disp('Unrecognized OS')
    return
end
filename = strcat('Li', mion, sprintf('F4_MC_%ux%ux%u_hyp=%u_', dims(1), dims(2), dims(3), hyperfineOpt), iState, '.mat');
params.DataFile = fullfile(filepath, char(filename));
params.saveOpt = saveOpt;
while true
    if loadOpt % continue simulation from previous file
        if isfile(params.DataFile)
            load(params.DataFile, "-mat", 'const', 'init', 'ion', 'params', 'wav', 'meas');
            const.gfac = const.mu0/4/pi * (ion.gLande(ion.idx) * const.muB)^2 * const.J2meV; % dipole interaction prefactor
            const.Ngfac = -ion.nLande(ion.idx) * const.muN * const.J2meV; % Nuclear Zeeman interaction prefactor
            coef = wav;
            if exist('meas', 'var')
                eSpin0 = meas.eSpin; % intermediate electronic spin configuration after thermalization
                nSpin0 = meas.nSpin; % intermediate nuclear spin configuration
                E_si = meas.Esi; % initial single-ion energy
            else
                eSpin0 = init.eSpinT; % use the initial condition
                nSpin0 = init.nSpinT;
                E_si = init.Esi; % initial single-ion energy
            end
            break;
        else
            loadOpt = false;
            fprintf('Load file not found, starting a new simulation!\n');
        end
    else % start a simulation from scratch
        %% define sample

        % define constants
        const.muB = 9.274e-24; % [J/T]
        const.kB = 1.3806e-23; % [J/K]
        const.muN = 5.05078e-27; % [J/T]
        const.mu0 = 4e-7 * pi; % [H/m]
        const.J2meV = 6.24151e+21; % [mev/J]

        % Ions' names
        ion.name = [{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
        ion.name_hyp = [{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];
        ion.prop = [0; 0; 0; 0; 0; 0]; % ion proportion
        ion.cfRot = [0; -11; 0; 0; 0; 0]; % crystal field basis rotation (in degrees)
        params.hyp = hyperfineOpt; % hyperfine interaction option
        ion.alt = [1  1  1
            1  1  1
            1  1  1
            1  1  1]; % default ferromagnetic order parameter counting
        % Ions' proportions
        switch mion
            case 'Er'
                ion.prop(1) = 1; % LiErF4
                ion.idx = 1; % index of the Re ion for calculation
                ion.oParams = 'x';
                ion.alt = [1  1  1
                    -1  1  1
                    -1 -1  1
                    1 -1  1]; % alternate counting order parameter
            case 'Ho'
                ion.prop(2) = 1; % LiHoF4
                ion.idx = 2;
                ion.oParams = 'z';
            case 'Yb'
                ion.prop(3) = 1; % LiYbF4
                ion.idx = 3;
                ion.oParams = 'x';
            case 'Tm'
                ion.prop(4) = 1; % LiTmF4
                ion.idx = 4;
                ion.oParams = 'z';
            case 'Gd'
                ion.prop(5) = 1; % LiGdF4
                ion.idx = 5;
            case 'Y'
                ion.prop(6) = 1; % LiYF4
                ion.idx = 6;
            case 'dope' % doped compound
                ion.prop(6) = 0.3; % doping percentage [Er, Ho, Yb, Tm, Gd, Y]
                ion.prop(2) = 1-ion.prop(1);
                ion.idx = [2 6];
        end

        % Concentrations of isotopes with nuclear spins
        ion.hyp = [0.23; 1; 0; 0; 0.14; 0]; % Er, Ho, Yb, Tm, Gd, Y
        % ion.hyp = [1; 1; 0; 0; 0; 0]; % Er, Ho, Yb, Tm, Gd, Y
        if ion.hyp(ion.idx) == 0; params.hyp = 0; end % Switch off hyperfine interaction depending on the magnetic ion

        % Ions' J, L, S, I values
        %        Er       Ho      Yb       Tm       Gd      Y
        ion.J = [15/2;    8;      7/2;    6;      7/2;    1];
        ion.L = [6;       6;      3;      5;      0;      1];
        ion.S = [3/2;     2;      1/2;    1;      7/2;    1];
        ion.I = [3;      3.5;      0;     0;      1.50;   0];
        ion.nLande = [-0.1611; 1.192; -0.2592; -0.462; -0.2265; 0]; % https://easyspin.org/documentation/isotopetable.html
        ion.gLande = gLande(ion.L, ion.S);

        % Hyperfine coupling strength
        A_Ho = 3.361e-3; % Ho hyperfine coupling (meV) -- Phys. Rev. B 75, 054426 (2007)
        % A_Ho = 1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000; % 1985 paper
        A_Er = 0.00043412; % Er hyperfine coupling (meV) -- Phys. Rev. B 2, 2298 - 2301 (1970)
        A_Gd = A_Ho;
        ion.A = [A_Er; A_Ho; 0; 0; A_Gd; 0];

        % Ions' lattice parameters [Angstrom]
        ion.abc = [{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Er
            {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      % Ho
            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      % Yb
            {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      % Tm
            {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Gd
            %            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}     % Gd
            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    % Y
        % basis vector
        ion.tau = [ 0     0     0
            0    1/2   1/4
            1/2   1/2   1/2
            1/2    0    3/4];

        % Ions' cf parameters (ueV)
        ion.B = [[60.2258   -0.1164   -4.3280   0.00   -0.0019   -0.0850   -0.0227]   % Er
            [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]          % Ho -- Phys. Rev. B 75, 054426 (2007)
            % [-57.9   0.309   3.51   0.00   0.000540   0.0631   0.0171]          % Ho -- Phys. Rev. B 92, 144422 (2015)
            [646.2016   15.3409  116.4854   0.00   -0.0686  -15.1817    0.0000]  % Yb
            %         [663.0000   12.5000   102.0000   0.00   -0.6200   -16.0000   0.0000] % Yb
            [224.3   -1.85   -11.7   0.00   0.002   0.2645   0.1377]             % Tm
            [0 0 0 0 0 0 0]                                                 % Gd
            [0 0 0 0 0 0 0]];                                               % Y
        ion.B = ion.B/1000.0; % ueV to meV

        ion.Hcf  = cf(ion.J(ion.idx), ion.B(ion.idx,:), ion.cfRot(ion.idx)); % CEF hamiltonian
        % O44 = (Jp^4 + Jm^4)/2; % off-diagonal CEF term for in-plane anisotropy
        % H_h4 = ion.h4(ion.idx) * O44 * h_int(1)^2; % order by disorder anisotropy

        % phenomenological anisotropy term
        % ion.h4 = [0.11e-3; 0; 0; 0; 0; 0]; % Crystal field anisotropy, Neda's thesis
        ion.h4 = [5.26e-3; 0; 0; 0; 0; 0]; % Order-by-disorder anisotropy, Neda's thesis

        % Exchange interaction strength
        ex_Er = 0;
        ex_Ho = 1e-4; % [meV] AFM exchange [PRB 75, 054426 (2007)]
        ex_Yb = 0;
        ex_Tm = 0;
        ion.ex = [ex_Er; ex_Ho; ex_Yb; ex_Tm; 0; 0];

        J = ion.J(ion.idx);
        ion.Jz = diag(J:-1:-J);
        Jp = diag(sqrt((J - (J-1:-1:-J)).*(J+1 + (J-1:-1:-J))),1);
        Jm = Jp';
        ion.Jx = (Jp+Jm)/2;
        ion.Jy = (Jp-Jm)/2i;

        % initialize the nuclear spins
        if params.hyp
            I = ion.I(ion.idx);
            ion.Iz = diag(I:-1:-I);
            Ip = diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1);
            Im = Ip';
            ion.Ix = (Ip + Im)/2;
            ion.Iy = (Ip - Im)/2i;
        end
        const.gfac = const.mu0/4/pi * (ion.gLande(ion.idx) * const.muB)^2 * const.J2meV; % dipole interaction prefactor
        const.Ngfac = -ion.nLande(ion.idx) * const.muN * const.J2meV; % Nuclear Zeeman interaction prefactor

        %% Calculation parameters
        switch mion
            case 'Er'
                % For LiErF4
                Bfield = 0;
                params.temp(1:5) = linspace(0.02,0.06,5);
                params.temp(6:27) = linspace(0.07,0.16,22);
                params.temp(28:32) = linspace(0.18,0.25,5);
            case 'Yb'
                % For LiYbF4
                Bfield = 0;
                params.temp(1:12) = linspace(0.01,0.12,12);
                params.temp(13:27) = linspace(0.120,0.140,15);
                params.temp(28:32) = linspace(0.142,0.15,5);
            case 'Ho'
                % For LiHoF4
                params.temp = 0.01;
                Bfield = linspace(0,6,24);
            otherwise
                prompt = sprintf('Please specify temperature range:\n');
                params.temp = input(prompt);
                prompt = sprintf('Please specify magnetic field range:\n');
                Bfield = input(prompt);
        end

        % % Log spaced temperature points
        % params.temp = zeros(1,40);
        % logtemp = logspace(0,-5,25);
        % params.temp(1:25) = 0.130*(1-logtemp);
        % params.temp(26:40) = linspace(0.130,0.250,15);

        % % Linear spaced temperature points
        % params.temp = zeros(1,40);
        % params.temp(1:7) = linspace(0,0.016,7);
        % params.temp(8:33) = linspace(0.020,0.040,26);
        % params.temp(34:40) = linspace(0.044,0.080,7);

        % magnetic field orientation
        theta = 0; % out-of-plane angle (radian) from ab plane
        phi = 0; % in-plane angle (radian) from a/b axis
        Bx = Bfield * cos(phi) * cos(theta);
        By = Bfield * sin(phi) * cos(theta);
        Bz = Bfield * sin(theta);
        params.field = [Bx; By; Bz];
        params.theta = theta;
        params.phi = phi;

        % initialization option
        params.init = iState; % initial spin configuration

        % % sample size
        % params.convg = 1e-9; % Energy convergence criteria for thermalization stage
        % params.tEQ = 1e4; % Thermalization steps (in unit of lattice size)
        % params.sampSize = 50; % Sampling size (in unit of lattice size)
        % params.mIntv = 1e3; % Interval between measurements
        % params.sIntv = 10; % data saving interval
        % % params.pt_intv = 100; % Interval between parallel temperature trials

        % sample size for debugging
        params.convg = 1e-9; % Energy convergence criteria for thermalization stage
        params.tEQ = 1e3; % Thermalization steps (in unit of lattice size)
        params.sampSize = 5; % Sampling size (in unit of lattice size)
        params.mIntv = 2; % Interval between measurements
        params.sIntv = 3; % Data saving interval
        % params.pt_intv = 1; % Interval between parallel temperature trials

        % construct lattice
        params.dims = dims; % N x N x N unit cells
        params.coord = 'spherical'; % sample shape (spherical, cubic, ellipsoidal)
        params.domRng = max(dims);  % domain cutoff (in unit of unit cells)
        params.dpRng = inf;  % domain cutoff & dipole summation range (in unit of unit cells)
        params.pos = lattice(ion.abc{ion.idx}, ion.tau, params); % lattice site coordinates

        params.temp = sort(params.temp); % order the temperature points from low to high
        switch iState
            case 'disordered'
                params.temp = flip(params.temp); % start at high temperature and anneal the ensemble
                coords = randSph(size(params.pos,1),'Cartesian');
                [alp, bet, ~] = cart2sph(coords(:,1),coords(:,2),coords(:,3));
                alp = alp + pi; % [-pi pi] -> [0 2pi]
                bet = bet + pi/2; % [-pi/2 pi/2] -> [0 pi]
                coef = [cos(bet'/2) ; sin(bet'/2).*exp(1i*alp')];
            case 'ordered'
                switch mion
                    case {'Ho', 'Y', 'Gd', 'Tm'} % z-FM ordered state
                        alp = zeros(1,size(params.pos,1));
                        bet = zeros(1,size(params.pos,1));
                        coef = [cos(bet/2); sin(bet/2).*exp(1i*alp)];
                    case {'Er', 'Yb'} % xy-AFM ordered state
                        alp = zeros(1,size(params.pos,1));
                        % alp = ones(1,size(params.pos,1)) * pi/2;
                        bet = ones(1,size(params.pos,1)) * pi/2;
                        coef = [cos(bet/2); sin(bet/2).*exp(1i*alp)];

                        % use the propogation factor to order the spins
                        k_afm = [1 0 0]' .* (2*pi ./ diag(ion.abc{ion.idx})); % X propogation
                        % k_afm = [0 1 0]' .* (2*pi ./ diag(ion.abc{ion.idx})); % Y propogation
                        coef(2,:) = coef(2,:) .* exp(1i * (params.pos * k_afm))'; % Apply the propogation vector
                end
            otherwise
                disp('Unrecognized initial condition!\n')
        end
        coef = repmat(coef, 1, 1, length(params.temp), size(params.field,2));
        break;
    end
end

% project the electronic hamiltonian to Ising space
[~, hamI, basis, ~, ~, ~] = Ising_proj(const, ion, params); % Ising projection
hamI = gpuArray(hamI);
coef = gpuArray(coef);
basis = gpuArray(basis);

if  ~loadOpt
    % Initialize the configuration for new simulations
    [E_si, eSpin0, nSpin0, params] = initialize_GPU(const, ion, params, coef, basis, hamI);
    if saveOpt == true
        save(params.DataFile, 'const', 'params', 'ion', '-v7.3');
    end
end

[init, coef, E_si] = equilibrate_GPU(const, ion, params, E_si, hamI, coef, basis, eSpin0, nSpin0);
wav = gather(coef);
init.eSpin0 = gather(eSpin0); % initial electronic spin configuration
init.nSpin0 = gather(nSpin0);% initial nuclear spin configuration
if saveOpt == true
    save(params.DataFile, 'wav', 'init', '-append');
end

% Monte Carlo sampling
[coef, meas] = MC_sample_GPU(const, ion, params, E_si, hamI, basis, wav, init);

if params.saveOpt == true
    save(params.DataFile, 'coef', 'meas', '-append');
    fprintf('Data saved!\n');
end

if plotOpt
    if ~loadOpt
        plot_spin(ion, init, meas, params, [0 0], {'thermalization'}, 'electron')
    end
    if length(params.temp) > size(params.field,2)
        plot_spin(ion, init, meas, params, [0 0], {'T_mag', 'vector', [ion.oParams, 'domain']}, 'electron')
    else
        plot_spin(ion, init, meas, params, [0 0], {'B_mag', 'vector', [ion.oParams, 'domain']}, 'electron')
    end
end
end