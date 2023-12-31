function [init, meas, coef, ion, params] = LiReF4_CMC_Yikai_a(mion, dim, hyperfineOpt, loadOpt, plotOpt, saveOpt)
% Argument: mion, loadopt, dim, loadOpt, plotOpt, saveOpt
% mion: Magnetic ion type
% dim: lattice dimentions (1 x 3 array)
% loadOpt: initiation option: 1. load existing file, 2. start from scratch
% plotOpt: plotting option: True/False
% saveOpt: saving option: True/False

% define constants
const.muB = 9.274e-24; % [J/T]
const.kB = 1.3806e-23; % [J/K]
const.muN = 5.05078e-27; % [J/T]
const.mu0 = 4e-7 * pi; % [H/m]
const.J2meV = 6.24151e+21; % [mev/J]

% set base folder depending on the OS
if ispc
    filepath = ['G:\My Drive\File sharing\PhD program\Research projects\Li', mion,...
        'F4 project\Data\Simulations\MATLAB\Monte Carlo'];
elseif ismac
    filepath = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
        'PhD program/Research projects/Li', mion, 'F4 project/Data/Simulations/Matlab/Monte Carlo/'];
else
    disp('Unrecognized OS')
    return
end
filename = strcat('Li', mion, sprintf('F4_MC_%ux%ux%u_hyp=%u', dim(1), dim(2), dim(3), hyperfineOpt), '.mat');
params.DataFile = fullfile(filepath, char(filename));

while true
    switch loadOpt
        case 'new' % start a simulation from scratch
            %% define sample
            % Ions' names
            ion.name = [{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
            ion.name_hyp = [{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];
            ion.prop = [0; 0; 0; 0; 0; 0]; % ion proportion
            ion.cfRot = [0; -11; 0; 0; 0; 0]; % crystal field basis rotation (in degrees)
            params.hyp = hyperfineOpt; % hyperfine interaction option

            % Ions' proportions
            switch mion
                case 'Er'
                    ion.prop(1) = 1; % LiErF4
                    ion.idx = 1; % index of the Re ion for calculation
                case 'Ho'
                    ion.prop(2) = 1; % LiHoF4
                    ion.idx = 2;
                case 'Yb'
                    ion.prop(3) = 1; % LiTbF4
                    ion.idx = 3;
                case 'Tm'
                    ion.prop(4) = 1; % LiTmF4
                    ion.idx = 4;
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
            % O44 = (Jp^4+Jm^4)/2; % off-diagonal CEF term for in-plane anisotropy
            % H_h4 = ion.h4(ion.idx) * O44 * h_int(1)^2; % order by disorder anisotropy

            % phenomenological anisotropy term
            % ion.h4 = [0.11e-3; 0; 0; 0; 0; 0]; % Crystal field anisotropy, Neda's thesis
            ion.h4 = [5.26e-3; 0; 0; 0; 0; 0]; % Order-by-disorder anisotropy, Neda's thesis

            % Heisenberg exchange interaction
            ex_Er = 0;
            ex_Ho = -1e-4; % [meV] AFM exchange [notation convention: PRB 75, 054426 (2007)]
            ex_Yb = 0;
            ex_Tm = 0;
            ion.ex = [ex_Er; ex_Ho; ex_Yb; ex_Tm; 0; 0];

            J = ion.J(ion.idx);
            ion.Jz = diag(J:-1:-J);
            Jp = diag(sqrt((J - (J-1:-1:-J)).*(J+1 + (J-1:-1:-J))),1);
            Jm = Jp';
            ion.Jx = (Jp+Jm)/2;
            ion.Jy = (Jp-Jm)/2i;

            %% Calculation setting
            % params.temp(1:6) = linspace(0.1,0.6,6); % for debugging
            % params.temp(1:4) = linspace(0.05,0.1,4);
            % params.temp(5:28) = linspace(0.11,0.15,24);
            % params.temp(29:32) = linspace(0.16,0.25,4);
            params.temp = 0.01;
            % params.temp = zeros(1,40);
            % % logtemp = logspace(0,-5,25);
            % % params.temp(1:25) = 0.130*(1-logtemp);
            % params.temp(26:40) = linspace(0.130,0.250,15);
            % params.temp = zeros(1,40);
            % params.temp(1:7) = linspace(0,0.016,7);
            % params.temp(8:33) = linspace(0.020,0.040,26);
            % params.temp(34:40) = linspace(0.044,0.080,7);

            % magnetic field (T)
            % Bfield = [0 5 10];
            Bfield = linspace(0,6,24);
            % Bfield(1:15) = linspace(0,0.8,15)';
            % Bfield(16:25) = linspace(0.84,19.8,10)';
            % Bfield(26:40) = linspace(20,24,15)';

            % magnetic field orientation
            theta = 0; % convert to radian, out-of-plane angle from ab plane
            phi = 0; % convert to radian, in-plane angle from a/b axis
            Bx = Bfield * cos(phi) * cos(theta);
            By = Bfield * sin(phi) * cos(theta);
            Bz = Bfield * sin(theta);
            params.field = [Bx; By; Bz];
            params.theta = theta;
            params.phi = phi;

            % initialization option
            params.init = 'random'; % initial spin configuration

            % sample size
            params.convg = 1e-7; % Energy convergence criteria for thermalization stage
            params.tEQ = 1e3; % Thermalization steps (in unit of lattice size)
            params.tMeas = 50; % Sampling size (in unit of lattice size)
            params.mIntv = 1e3; % Interval between measurements
            params.sIntv = 10; % data saving interval
            % params.pt_intv = 100; % Interval between parallel temperature trials

%             % for debugging
%             params.convg = 1e-7; % Energy convergence criteria for thermalization stage
%             params.tEQ = 1e3; % Thermalization steps (in unit of lattice size)
%             params.tMeas = 10; % Sampling size (in unit of lattice size)
%             params.mIntv = 100; % Interval between measurements
%             params.sIntv = 5; % data saving interval
%             % params.pt_intv = 1; % Interval between parallel temperature trials

            % construct lattice
            params.dims = dim; % N x N x N unit cells
            params.coord = 'spherical'; % sample shape (spherical, cubic, ellipsoidal)
            params.domRng = max(dim);  % domain cutoff & dipole summation range (in unit of unit cells)
            params.dpRng = 2 * sqrt(max(params.dims)^3);  % domain cutoff & dipole summation range (in unit of unit cells)
            params.pos = lattice(ion.abc{ion.idx}, ion.tau, params); % lattice site coordinates

            if strcmpi(params.init, 'random')
                coords = randSph(size(params.pos,1),'Cartesian');
                [alp, bet, ~] = cart2sph(coords(:,1),coords(:,2),coords(:,3));
                alp = alp + pi; % [-pi pi] -> [0 2pi]
                bet = bet + pi/2; % [-pi/2 pi/2] -> [0 pi]
                coef = [cos(bet'/2) ; sin(bet'/2).*exp(1i*alp')];
            elseif strcmpi(params.init, 'uniform')
                alp = zeros(1,size(params.pos,1));
                bet = zeros(1,size(params.pos,1));
                coef = [cos(bet/2); sin(bet/2).*exp(1i*alp)];
            else
                disp('Unrecognized initial condition!\n')
            end
            coef = repmat(coef, 1, 1, length(params.temp), size(params.field,2));
            break;
        case 'load' % continue simulation from previous file
            if isfile(params.DataFile)
                load(params.DataFile, "-mat", 'init', 'ion', 'params', 'coef', 'meas');
                eSpinT = meas.eSpin; % intermediate spin configuration after thermalization
                Esi = meas.en0; % initial single-ion energy
                break;
            else
                loadOpt = 'new';
                fprintf('Load file not found, starting a new simulation!\n');
            end
        otherwise
            fprintf('Loading option unspecified, starting a new simulation!\n');
            loadOpt = 'new';
    end
end

% project the electronic hamiltonian to Ising space
[~, hamI, basis, ~, ~, ~] = Ising_proj(const, params.temp, params.field, ion); % Ising projection

% Initialize the configuration for new simulations
if  strcmp(loadOpt, 'new')
    eSpin0 = double.empty(size(params.pos,1), 3, length(params.temp), size(params.field, 2), 0); % initial electornic spin config
    Esi = double.empty(size(params.pos,1), length(params.temp), size(params.field, 2), 0); % single-ion energy
    Edip = double.empty(size(params.pos,1), length(params.temp), size(params.field, 2), 0); % global energy update trace
    nSpin0 = zeros(size(params.pos,1), 3, length(params.temp), size(params.field, 2)); % final spin config
    Enuc = zeros(size(params.pos,1), length(params.temp), size(params.field, 2)); % global energy update trace
    for tt = 1:length(params.temp)
        for ff = 1:size(params.field,2)
            bas = squeeze(basis(:,:,tt,ff)); % single-ion wavefunction
            for ii = 1:size(params.pos,1)
                % initialize the electronic spins
                wav = squeeze(coef(:,ii,tt,ff)); % Unitary rotation
                spx = real(wav' * bas' * ion.Jx * bas * wav)';
                spy = real(wav' * bas' * ion.Jy * bas * wav)';
                spz = real(wav' * bas' * ion.Jz * bas * wav)';
                eSpin0(ii,:,tt,ff,1) = [spx, spy, spz]; % initial electronic spin config
                Esi(ii,tt,ff,1) = real(wav' * squeeze(hamI(:,:,tt,ff)) * wav); % initial single-ion energy
               
                if params.hyp
                    % initialize the nuclear spins
                    I = ion.I(ion.idx);
                    ion.Iz = diag(I:-1:-I);
                    Ip = diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1);
                    Im = Ip';
                    ion.Ix = (Ip + Im)/2;
                    ion.Iy = (Ip - Im)/2i;

                    H_hyp = ion.A(ion.idx) * (spx*ion.Ix + spy*ion.Iy + spz*ion.Iz); % nuclear hyperfine interaction (meV)
                    HzI = -ion.nLande(ion.idx) * const.muN * const.J2meV * (Bfield(1)*ion.Ix +...,
                        Bfield(2)*ion.Iy + Bfield(3)*ion.Iz); % nuclear Zeeman interaction (meV)
                    ham_nuc = H_hyp + HzI; % total nuclear spin hamiltonian
                    [eigen_n, ~] = eig(ham_nuc); % diagonalize the nuclear hamiltonian
                    wav_n = randSphN(1,size(ham_nuc,1)); % random point of the nuclear spin configuration
                    spx_n = real(wav_n' * eigen_n' * ion.Ix * eigen_n * wav_n)';
                    spy_n = real(wav_n' * eigen_n' * ion.Iy * eigen_n * wav_n)';
                    spz_n = real(wav_n' * eigen_n' * ion.Iz * eigen_n * wav_n)';
                    nSpin0(ii,:,tt,ff,1) = [spx_n spy_n spz_n]; % initial nuclear spin config
                    Enuc(ii, tt,ff,1) = real(wav_n' * eigen_n' * ham_nuc * eigen_n * wav_n); % new energy
                    Esi(ii,tt,ff,1) = Esi(ii,tt,ff,1) + Enuc(ii, tt,ff,1); % include hyperfine and nuclear Zeeman energies
                end
            end
            % calculate the initial dipole interaction energies
            for ii = 1:size(params.pos,1)
                Edip(ii, tt,ff, 1) = CntrDip(const, ion, params.pos, eSpin0(:,:,tt,ff), ii, eSpin0(ii,:,tt,ff)); % initial dipolar interaction energy
            end
        end
    end
    fprintf('Initialization compelte!\n')

    [eSpinT, nSpinT, coef, Esi, Edip, Enuc, dE, accpRate] = equilibrate_a(const, ion, params, Esi, Edip, hamI, coef, basis, eSpin0, Enuc, nSpin0);
    init.eSpin0 = eSpin0; % initial electronic spin configuration
    init.eSpinT = eSpinT; % intermediate electronic spin configuration after thermalization
    init.nSpin0 = nSpin0;% initial nuclear spin configuration
    init.nSpinT = nSpinT; % intermediate nuclear spin configuration after thermalization
    init.dEt = dE; % history of energy change during thermalization
    init.aRate = accpRate; % history of acceptance rate during thermalization
    if saveOpt == true
        save(params.DataFile, 'init', 'params', 'ion', '-v7.3');
    end
end

% Monte Carlo sampling
SampSize = 0;
while SampSize < params.tMeas
    [Mx, My, Mz, coef, eSpinT, nSpinT, Esi, Edip, ~] = MC_sample_a(const, ion, params, Esi, Edip, hamI, basis, coef, eSpinT, nSpinT, Enuc);
        if SampSize == 0
            meas.Mx = Mx;
            meas.My = My;
            meas.Mz = Mz;
        else
            meas.Mx = cat(1, meas.Mx, Mx);
            meas.My = cat(1, meas.My, My);
            meas.Mz = cat(1, meas.Mz, Mz);
        end
        meas.eSpin = eSpinT;
        meas.nSpin = nSpinT;
        meas.en0 = Esi + Edip;
    if saveOpt == true
        save(params.DataFile, 'coef', 'meas', '-append');
        fprintf('Data saved!\n');
    end
    SampSize = SampSize + params.sIntv;
end

if plotOpt
    if strcmp(loadOpt, 'new')
        fg1 = figure;
        ax1 = fg1.CurrentAxes;
        lgd = [];
        hold on;
        for tt = 1:length(params.temp)
            for ff = 1:size(params.field, 2)
                % Plot with specific color/line style and store the handle
                dEt = squeeze(dE(:, tt, ff));
                cutoff = find(dEt,1,'last');
                pl(tt, ff) = plot(1:cutoff, dEt(1:cutoff), 'LineWidth', 1.5);
                lgd{end+1} = ['Temp: ' num2str(params.temp(tt)) ', Field: ' num2str(ff)];
            end
        end
        xlabel('Iteration', 'Interpreter', 'latex');
        ylabel('$\Delta E$', 'Interpreter', 'latex');
        title('Energy Change over Iterations', 'Interpreter', 'latex');
        set(ax1, 'FontSize', 12);
        grid on;
        legend(lgd, 'Location', 'best');
    end

    fg2 = figure;
    ax2 = fg2.CurrentAxes;
    plot(vecnorm(params.field,1), squeeze(rms(meas.Mx,1)), 'o');
    hold on
    plot(vecnorm(params.field,1), squeeze(rms(meas.My,1)), 'o');
    plot(vecnorm(params.field,1), squeeze(rms(meas.Mz,1)), 'o');
    legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
    xlabel('Magnetic Field (T)')
    ylabel('$\langle J \rangle$','Interpreter','latex')
    set(ax2, 'FontSize', 12);

    plot_spin(init, meas, params, 1, {'vector', 'domain'}, 'electron')
end

if saveOpt
    save(params.DataFile, 'coef', 'meas', '-append')
end
end