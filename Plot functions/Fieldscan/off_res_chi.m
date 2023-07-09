Options.RPA  = true;
Options.filfctr = 0.0114; % filling factor correction
% Options.filfctr = 0.112; % filling factor correction
Options.Bcorrt = 4.4/4.6; % field correction (default = 1)
Options.nZee = true; % Nuclear Zeeman interaction
Options.search = 'automatic'; % 1. 'Automatic' search and match or, 2. 'manual' appointment of files 
    sample = 239; % Sample to use (SC###)
    sim_temps = [0.15 0.3 0.5 0.8]; % temperatures
    % sim_temps = 0.15;
    theta = 0.0; % c-axis deviation
    phi = 15.0; % in-plane (a-b) rotation
    freq_l = 4.69; % frequency range lower limit [GHz]
    freq_h = 4.745; % frequency range upper limit [GHz]
    gama = 9e-5; % default spin linewidth
    if strcmp(pathsep, ':')
        platform = 'Unix';
    else
        platform = 'Win';
    end
switch Options.search
    case 'automatic'
        % set experimental data loading parameters
        switch platform
            case 'Win' % for Windows
                exp_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\Research projects\',...
                    'LiHoF4 project\Data\Experiment\LiHoF4\SC', num2str(sample), '\'];
                sim_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\',...
                    'Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities\'];
                if Options.nZee == true
                    sim_loc = [sim_loc,'Hz_I=1\'];
                else
                    sim_loc = [sim_loc, 'Hz_I=0\'];
                end
            case 'Unix' % for Mac/Linux
                exp_loc = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
                    'PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC', num2str(sample), '/'];
                sim_loc = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
                    'PhD program/Research projects/LiHoF4 project/Data/Simulations/MATLAB/Susceptibilities/'];
                if Options.nZee == true
                    sim_loc = [sim_loc,'Hz_I=1/'];
                else
                    sim_loc = [sim_loc, 'Hz_I=0/'];
                end
        end
        
        dat_lst = readtable(fullfile(exp_loc, sprintf("SC%u_off_list.xlsx", sample)));
        folder = table2array(dat_lst(:,1));
        efiles = table2array(dat_lst(:,2));
        etemps = table2array(dat_lst(:,3));
    case 'manual'
        switch platform
            case 'Win'
                folder = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\',...
                    'File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
                    'SC107\19.05.2019\']; % for Windows
            case 'Unix'
                folder = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
                    'PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/',...
                    'SC107/19.05.2019/']; % for Mac/Linux
                efiles = '2019_05_0026';
        end
end

muB = 9.274e-24; % Bohr magneton [J/T]
mu0 = 4*pi*1e-7; % Vacuum permeability ([H/m])
rho = 4/(5.175e-10 * 5.175e-10 * 10.750e-10); % magnetic moment number density in LiReF4 [m^-3]
hbar = 1.05457E-34; % Reduced Planck constant [J.s/rad]
J2meV = 6.24151e+21; % [meV/J]
Gh2mV = hbar * 2*pi * 10^9 * J2meV; % [meV/GHz]
ELEf = gLande(6,2) * muB * J2meV; % [meV/T]

mkr = ['o','s','^','x','<','>']; % marker style pool
lin = ["-", "--" , "-.", ":"]; % line style pool
mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool

lgd_sim = string(sim_temps.*1000); % legends for calculated data
lgd_sim = 'Sim.' + lgd_sim + ' mK';
lgd_sim = lgd_sim';
% lgd_sim = strings(length(temps),1); % container for legends for simulated data
exp_temps = double.empty(length(sim_temps),0); % container for legends for experimental data

% plot the data
chir_fig = figure;
box on
hold on
chii_fig = figure;
box on
hold on
for ii = 1:length(sim_temps)
    switch Options.search
        case 'automatic'
            [~, idx] = min(abs(etemps-sim_temps(ii))); % find the closest experimental temperature
            exp_file = fullfile([exp_loc, folder{idx}], efiles{idx}); % locate the experimental data
            exp_temps(ii,1) = etemps(idx); % record the actual experimental temperatures
            load([exp_file, '_interp'], '-mat', 'analysis', 'continu_var');
        case 'manual'
            exp_file = fullfile(folder, efiles);
            load([exp_file, '_interp'], '-mat', 'analysis', 'continu_var');
            exp_temps(ii,1) = analysis.temp;
    end
    f0 = analysis.w0; % resonant frequency
    B0 = analysis.Hs;
    wc = nonzeros(analysis.wc); % cavity resonant frequency [GHz]
    gw0 = sqrt(mu0 * 2*pi * wc(1)*1e9 * rho/2) * Options.filfctr; % susceptibility prefactor [T^2/J. rad/s]^1/2
    gw2 = gw0^2 * 2*pi * 1e-9 / J2meV; % [GHz * T^2/meV]
    
    if ~isfield(Options,'Bcorrt'); Options.Bcorrt = 1; end
    figure(chir_fig)
    scatter(continu_var(1,1:15:end)*Options.Bcorrt, analysis.xr(1:15:end)/gw2, 'Marker', mkr(mod(ii,length(mkr))),...
        'MarkerEdgeColor', mc{mod(ii,size(mc,1))}, 'MarkerFaceColor', mc{mod(ii,size(mc,1))}); % plot the experimental Re[chi]

    figure(chii_fig)
    scatter(continu_var(1,:)*Options.Bcorrt, medfilt1(analysis.xi)/gw2, 'Marker', '.',...
        'MarkerEdgeColor', mc{mod(ii,size(mc,1))}, 'MarkerFaceColor', mc{mod(ii,size(mc,1))}); % plot the experimental Im[chi]
%     figure(chii_fig)
%     plot(continu_var(1,1:15:end), analysis.xi(1:15:end), 'LineStyle', '-', 'Marker', mkr(mod(ii,length(mkr))),...
%         'Color', mc{mod(ii,size(mc,1))}); % plot the experimental data
    if exist('freq_l','var') && exist('freq_h','var')
        sim_file = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_%5$.3f-%6$.3fGHz',...
            sim_temps(ii), theta, phi, gama, freq_l, freq_h);
    else
        sim_file = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_*GHz', sim_temps(ii), theta, phi, gama);
    end
    sim_load = strcat("chi_LiHoF4_", sim_file, ".mat");
    chi_file = dir(fullfile(sim_loc, sim_load)); % search for the specified file
    if isempty(chi_file)
        disp(strcat("Cannot match the file: ", sim_load, ", check the name and the directory!"))
        break
    end
    if Options.RPA == true
        xlabr = 'Re[\chi^{zz}_{RPA}] (meV)';
        xlabi = 'Im[\chi^{zz}_{RPA}] (meV)';
        load([sim_loc chi_file.name],'-mat','freq_total','fields','chiq','unit'); % load calculated suscpetibilities
        chiz = squeeze(chiq(3,3,:,:)); % extract zz tensor element
    else
        xlabr = 'Re[\chi^{zz}_{MF}] (meV)';
        xlabi = 'Im[\chi^{zz}_{MF}] (meV)';
        load([sim_loc chi_file.name],'-mat','freq_total','fields','chi','unit'); % load calculated suscpetibilities
    chiz = squeeze(chi(3,3,:,:)); % extract zz tensor element
    end

    % Convert the unit to meV
    if strcmp(unit, 'GHz')
        ConvUnit = Gh2mV;
    elseif strcmp(unit, 'J')
        ConvUnit = J2meV;
    elseif strcmp(unit, 'meV')
        ConvUnit = 1;
    else
        disp('Energy unit unrecognized!')
    end
    
    % use resonant frequency trace for linear plot of susceptibility
    fref = interp1(B0, f0, fields); % cavity resonance frequency
    % fref = 4.9164.*ones(size(fields)); % fixed cavity resonance frequency
    fidx = uint32.empty(length(fields),0); % resonant frequency indices
    chir = double.empty(length(fields),0);
    chii = double.empty(length(fields),0);
    for jj = 1:length(fields)
        [~, fidx(jj)] = min(abs(freq_total - fref(jj)));
        chir(jj) = ELEf^2 .* real( chiz(fidx(jj),jj) ) ./ ConvUnit;
        chii(jj) = ELEf^2 .* imag( chiz(fidx(jj),jj) ) ./ ConvUnit;
    end
%     [~, fidx] = min(abs(freq_total - fc)); % use simple cavity frequency
    
    figure(chir_fig);
    plot(fields, chir, 'LineStyle', lin(1), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);
%     plot(fields, real(chiz(fidx,:))/ConvUnit, 'LineStyle', lin(1), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);

    figure(chii_fig);
    plot(fields, chii, 'LineStyle', lin(2), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);
%     plot(fields, imag(chiz(fidx,:))/ConvUnit, 'LineStyle', lin(1), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);
end

lgd_exp = string(exp_temps.*1000); % legends for experimental data
lgd_exp = 'Exp.' + lgd_exp + ' mK';
lgd = [lgd_exp lgd_sim]';
lgd = lgd(:); % interweaved legends

figure(chir_fig);
legend(lgd, 'Location', 'northwest')
xlim([0 9])
title('Real part of AC susceptibilities (easy axis)')
xlabel('Magnetic Field (T)');
ylabel(xlabr);

figure(chii_fig);
legend(lgd, 'Location', 'northwest')
xlim([0 9])
title('Imaginary part of AC susceptibilities (easy axis)')
xlabel('Magnetic Field (T)');
ylabel(xlabi);
