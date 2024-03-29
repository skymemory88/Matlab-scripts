function [vvv, eee, fff, ttt] = LiReF4_MF_Yikai(mion,temp,Bfield,theta,phi,hyp)
% Argument: mion, temp, field, theta, phi
% mion: Magnetic ion type
% Temp(s): can be a single value or an array
% Field(s): can be a single value or an array
% theta(in degrees): angle between the field and a-axis in ac-plane, 0 means a perfect transverse magnetic field,
% phi(in degrees): ab-plane rotation, phi(in radian) = 0 means H along x

%% Setup parameters.
global rundipole muB mu0 muN J2meV
rundipole = true;
muB = 9.274e-24; % [J/T]
muN = 3.15245e-5; % [meV/T]
mu0 = 4e-7*pi; % [H/m]
J2meV = 6.24151e+21; % [mev/J]

global Options;
Options.hyperfine = hyp; % hyperfine interaction option
Options.nZee = true; % Nuclear Zeeman interaction option
Options.demag = true; % Demagnetization factor option
ion.alpha = 0; % shape in calculating demagnetization factor is a needle (0), sphere (1), disc (inf)
Options.plotting = false; % Choose whether or not to plot the figures at the end

global strategies; % global convergence strategies switches
strategies.powerlaw = false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator = 0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping = 0.2; % damping factor. May reduce exponential fit efficiency
strategies.expfit = true; % turn on fitting to an exponential.
strategies.expfit_period = 30; % period between exponential fit.
strategies.expfit_deltaN = 5; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.
strategies.symmetry = true; % Copy the state of one site to its crystal symmetry equivalent sites

if ispc
        Options.filepath = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing',...
            '\PhD program\Research projects\Li', mion, 'F4 project\Data\Simulations\Matlab\Susceptibilities\'];
elseif ismac
        Options.filepath = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
            'PhD program/Research projects/Li', mion, 'F4 project/Data/Simulations/Matlab/Susceptibilities/'];
end

if Options.nZee == true && Options.hyperfine == true
    Options.filepath = [Options.filepath,'Hz_I=1'];
else
    Options.filepath = [Options.filepath,'Hz_I=0'];
end

if length(Bfield) > length(temp)
    Options.scanMode = 'field';
else
    Options.scanMode = 'temp';
end

Options.theta = theta;
Options.phi = phi;

theta = theta*pi/180; % convert to radian, deviation angle from normal vector
phi = phi*pi/180; % convert to radian, in-plane angle from a/b axis

Hx = Bfield*cos(phi)*cos(theta);
Hy = Bfield*sin(phi)*cos(theta);
Hz = Bfield*sin(theta);
fields = [Hx; Hy; Hz];
%% define LiReF4
%Ions' names
ion.name = [{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
ion.name_hyp = [{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];
ion.prop = [0; 0; 0; 0; 0; 0]; % ion proportion
ion.cfRot = [0; -11; 0; 0; 0; 0]; % crystal field basis rotation (in degrees)
ion.h4 = 0; % phenomenological in-plane anisotropy term
        
%Ions' hyperfine
if Options.hyperfine
    ion.hyp = [0.23; 1; 0; 0; 0.14; 0]; % Er, Ho, Yb, Tm, Gd, Y
%     ion.hyp = [1; 1; 0; 0; 0; 0]; % Er, Ho, Yb, Tm, Gd, Y
else
    ion.hyp = [0; 0; 0; 0; 0; 0]; % turn off hyperfine interactions for all elements
end

%Ions' proportions
switch mion
    case 'Er'
        ion.prop(1) = 1; % LiErF4
        strategies.symmetry = false; % override symmetry reflection for Ising magnet
    case 'Ho'
        ion.prop(2) = 1; % LiHoF4
        strategies.symmetry = true;
    case 'Yb'
        ion.prop(3) = 1; % LiTbF4
        strategies.symmetry = true;
    case 'Tm'
        ion.prop(4) = 1; % LiTmF4
    case 'Gd'
        ion.prop(5) = 1; % LiGdF4
    case 'Y'
        ion.prop(6) = 1; % LiYF4
    case 'dope' % doped compound
        ion.prop(1) = 0.3; % doping percentage [Er, Ho, Yb, Tm, Gd, Y]
        ion.prop(2) = 1-ion.prop(1);
end

%Ions' J, L, S, I values
%      Er       Ho      Yb      Tm      Gd      Y
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
ion.dpRng = 100; % dipole summation range (in unit of unit cells)

% Ions' lattice parameters
ion.abc = [{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Er
           {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      % Ho
           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      % Yb
           {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      % Tm
           {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      % Gd
%            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}     % Gd
           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    % Y
% basis vector
ion.tau = [ 0    0    0
            0   1/2  1/4
           1/2  1/2  1/2
           1/2   0   3/4];
       
% Ions' cf parameters (ueV)
ion.B = [[60.2258   -0.1164   -4.3280   0.00   -0.0019   -0.0850   -0.0227]   % Er
%          [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]          % Ho -- Phys. Rev. B 75, 054426 (2007)
         [-60.0   0.350   3.60   0.00   0.000400   0.0655   0.0098]          % Ho -- SC239 (Jz/1.3004), {phi=7, nZ=0} or {phi=12, nZ=1}, demag=0
%          [-57.9   0.309   3.51   0.00   0.000540   0.0631   0.0171]          % Ho -- Phys. Rev. B 92, 144422 (2015)
%          [-57.9   0.309   3.60   0.00   0.000540   0.0570   0.0130]          % Ho -- PRB 2015 mod(Jz*0.89)
%          [-56.2   0.325   3.61   0.00   0.000181   0.0758   0.0000]          % Ho -- Handbook phys. & chem. Rare-Earths V.23 (1996)
%          [-73.6   0.478   4.69   0.00   0.000100   0.0861   0.0118]          % Ho -- J. Magn. Magn. Magn. 15, 31 (1980)
%          [-52.2   0.323   3.59   0.00   0.000522   0.0685   0.0000]          % Ho -- Phys. Rev. B 19, 6564 (1979)
%          [-52.0   0.281   3.70   0.00   0.000700   0.0704   0.0000]          % Ho -- Opt. Spec-trosc. 44, 68 (1978)
%          [-65.0   0.426   4.53   0.00   0.000100   0.0855   0.0169]          % Ho -- Phys. Rev. B 12, 5315 (1975)
%          [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0060]          % Ho -- original code
         [646.2016   15.3409  116.4854   0.00   -0.0686  -15.1817    0.0000]  % Yb    
%         [663.0000   12.5000   102.0000   0.00   -0.6200   -16.0000   0.0000] % Yb
         [224.3   -1.85   -11.7   0.00   0.002   0.2645   0.1377]             % Tm
         [0 0 0 0 0 0 0]                                                 % Gd
         [0 0 0 0 0 0 0]];                                               % Y
ion.B = ion.B/1000.0; % ueV to meV

% phenomenological anisotropy term 
% ion.h4 = [0.11e-3; 0; 0; 0; 0; 0]; % Crystal field anisotropy, Neda's thesis
ion.h4 = [5.26e-3; 0; 0; 0; 0; 0]; % Order-by-disorder anisotropy, Neda's thesis

%Ions' renormalization factor
ion.renorm = [[1, 1, 1]         % Er
             [1, 1, 1/1.3004]     % Ho -- Phys. Rev. B 75, 054426 (2007), Phys. Rev. B 49, 11833 (1994)
%              [1, 1, 0.89]        % Ho -- For PRB 2015 mod
             [1, 1, 1]         % Yb
             [1, 1, 1]        % Tm
             [1, 1, 1]         % Gd
             [1, 1, 1]];       % Y

% Heisenberg exchange interaction
ex.Er = 0;
ex.Ho = -1e-4; % [meV] AFM exchange [notation convention: PRB 75, 054426 (2007)]
ex.Yb = 0;
ex.Tm = 0;
ion.ex = [ex.Er; ex.Ho; ex.Yb; ex.Tm; 0; 0];

%Ions' initial moments
% ion.mom(:,:,1) = [1 0 0; 1 0 0; -1 0 0; -1 0 0];        % Er (Original code--Yikai)
ion.mom(:,:,1) = [1 0 0; -1 0 0; -1 0 0; 1 0 0];       % Er Domain-1
% ion.mom(:,:,1) = [0 1 0; 0 1 0; 0 -1 0; 0 -1 0];       % Er Domain-2
ion.mom(:,:,2) = [0 0 1;  0 0 1;  0 0 1; 0 0 1];        % Ho
% ion.mom(:,:,2) = [3.473 -0.045 0.1;...
%                   3.473 -0.045 0.1;...
%                   3.473 -0.045 0.1;...
%                   3.473 -0.045 0.1];                   % Ho
ion.mom(:,:,3) = [1 0 0; -1 0 0; -1 0 0; 1 0 0];        % Yb
ion.mom(:,:,4) = [1 0 0; -1 0 0; -1 0 0; 1 0 0];        % Tm
ion.mom(:,:,5) = [1 0 0; 1 0 0; -1 0 0; -1 0 0];        % Gd
ion.mom(:,:,6) = [0 0 0; 0 0 0; 0 0 0; 0 0 0];          % Y
for ii = 1:length(ion.name)
%     ion.mom(:,:,ii) = ion.mom(:,:,ii) .* ion.J(ii); % put in real spin moment length
%     ion.mom(:,:,ii) = ion.mom(:,:,ii) .* ion.J(ii) .* repmat(ion.renorm(ii,:)',1,size(ion.mom,1))'; % put in real spin moment length
    ion.mom(:,:,ii) = ion.mom(:,:,ii) .* 5.51;
end
%% Calculate and save results inside the LiIonsF4 function
[ion,~,E,V] = LiIonsF4(ion,temp,fields,phi,theta);
%% Plot data
n = ion.prop~=0;
if Options.plotting == true
    switch Options.scanMode
        case 'field'
            J_H = zeros([3,size(fields,2),size(ion.name,1)]);
            Jnorm_H = zeros([size(fields,2),size(ion.name,1)]);
        case 'temp'
            J_H = zeros([3,size(temp,2),size(ion.name,1)]);
            Jnorm_H = zeros([size(temp,2),size(ion.name,1)]);
    end
    
    for ii=1:size(ion.name,1)
        J_H(:,:,ii) = squeeze(mean(ion.Js_hyp(:,:,:,ii),1)); % With hyperfine
%         J_H(:,:,i) = squeeze(mean(ion.Js(:,:,:,:,i),1)); % Without hyperfine
        Jnorm_H(:,ii) = squeeze(mean(ion.Jmom_hyp_norm(:,ii),1)); % with hyperfine
%         Jnorm_H(:,i) = squeeze(mean(ion.Jmom_norm(:,:,i),1)); % Without hyperfine
    end
    
    HH = zeros(1,size(fields,2)); % Calculate the norm of the applied magnetic field
    for ii = 1:size(fields,2)
        HH(ii) = norm(fields(:,ii));
    end
    
    %Plot spin moments
    switch Options.scanMode
        case 'field' % For field scan
            figure
            hold on
            plot(HH(1:10:end),Jnorm_H(1:10:end,n),'-o','Color', 'black', 'LineWidth', 2, 'MarkerSize',4);
            title(sprintf('Norm of spin moments for %s', char(ion.name(n))));
            xlabel('Magnetic field (|B| [T])');
            ylabel('Spin moment norm <J_n>');
            
            figure
            hold on
            plot(HH,J_H(:,:,n),'-','LineWidth', 2);
            plot(HH(1:10:end), J_H(:,1:10:end,n),'o','MarkerSize',4); % Add spaced markers to the <J> curves
            plot(HH(1:10:end), Jnorm_H(1:10:end,n),'-s','MarkerSize',4);
            title(sprintf('Expectation values of spin moments for %s', char(ion.name(n))));
            legend('<J_x>','<J_y>','<J_z>','<J_{norm}>');
            xlabel('Magnetic field (|B| [T])');
            ylabel('Spin moment <J>');
        case 'temp'
            figure
            hold on
            p1 = plot(temp,ion.altJmom(:,:,n),'-o','Color', 'black', 'LineWidth', 2, 'MarkerSize',4);
            title(sprintf('Alternating spin moments for %s', char(ion.name(n))));
            xlabel('Temperature (K)');
            ylabel('Alternating Spin moment <altJ_n>');
            
            figure
            hold on
            plot(temp,J_H(:,:,n),'-','LineWidth', 2);
            plot(temp,J_H(:,:,n),'o','MarkerSize',4); % Add spaced markers to the <J> curves
            plot(temp(1:10:end),Jnorm_H(1:10:end,n),'-s','MarkerSize',4);
            title(sprintf('Expectation values of spin moments for %s', char(ion.name(n))));
            legend('<J_x>','<J_y>','<J_z>','<J_{norm}>');
            xlabel('Magnetic field (|B| [T])');
            ylabel('Spin moment <J>');
    end
end
%% Save the data files
% if more than one temperatures provided, the data is saved insied LiIonF4() function
if length(temp) == 1 || size(Bfield,2) == 1 % for single external paramter sweep, save the data outside LiIonsF4()
    eee = squeeze(E);
    fff = fields;
    ttt = temp;
    vvv = squeeze(V);
    switch  Options.scanMode
        case 'field'
            elem = [ion.name(n)];
            tit = strcat('Hscan_','Li',elem(1:end),sprintf('F4_%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                temp, theta*180/pi, phi*180/pi,ion.hyp(n)));
        case 'temp'
            elem = [ion.name(n)];
            tit = strcat('Tscan_Li',elem(1:end),sprintf('F4_%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                Bfield, theta*180/pi, phi*180/pi,ion.hyp(n)));
    end
    fileobj = fullfile(Options.filepath,char(tit));
    save(fileobj,'ttt','fff','eee','vvv','ion','Options','-v7.3')
end