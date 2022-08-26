Options.RPA  = true;
Options.nZee = false; % Nuclear Zeeman interaction
Options.sample = 239; % Sample to use (SC###)

hbar = 1.05457E-34; % Reduced Planck constant [J.s/rad]
meV2J = 1.60218e-22; % [J/meV]
Gh2mV = hbar * 2*pi / meV2J * 10^9/1.05; % [meV/GHz]

temps = [0.1 0.3 0.5 0.8]; % temperatures
% temp = 1;
theta = 0.0; % c-axis deviation
phi = 2.3; % in-plane (a-b) rotation
freq_l = 4.700; % frequency range lower limit [GHz]
freq_h = 4.745; % frequency range upper limit [GHz]
gama = 1e-4; % default spin linewidth

mkr = ['o','s','^','x','<','>']; % marker style pool
lin = ['-','--','-.',':']; % line style pool
mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool

lgd_sim = string(temps.*1000); % legends for calculated data
lgd_sim = 'Sim.' + lgd_sim + ' mK';
lgd_sim = lgd_sim';
% lgd_sim = strings(length(temps),1); % container for legends for simulated data
exp_temps = double.empty(length(temps),0); % container for legends for experimental data

% set experimental data loading parameters
exp_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC',...
    num2str(Options.sample),'\'];

dat_lst = readtable(fullfile(exp_loc, sprintf("SC%u_off_list.xlsx", Options.sample)));
folder = table2array(dat_lst(:,1));
efiles = table2array(dat_lst(:,2));
etemps = table2array(dat_lst(:,3));

% set calculated data loading parameters
if Options.nZee == true
    sim_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\',...
        'Susceptibilities\with Hz_I\'];
else
    sim_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\',...
        'Susceptibilities\without Hz_I\CEF_2015_mod_R=0.785_phi=2.3\'];
end

% plot the data
chir_fig = figure;
box on
hold on
chii_fig = figure;
box on
hold on
for ii = 1:length(temps)
    [~, idx] = min(abs(etemps-temps(ii))); % find the closest experimental temperature
    exp_file = fullfile([exp_loc, folder{idx},'\'], efiles{idx}); % locate the experimental data
    load([exp_file, '_interp'], '-mat', 'analysis', 'continu_var');
    fc = nonzeros(analysis.wc);
    fc = fc(1); % pick zero-field resonant frequency as the reference frequency
    exp_temps(ii,1) = etemps(idx); % record the actual experimental temperatures


    figure(chir_fig)
    scatter(continu_var(1,1:15:end), analysis.xr(1:15:end), 'Marker', mkr(mod(ii,length(mkr))),...
        'MarkerEdgeColor', mc{mod(ii,size(mc,1))}, 'MarkerFaceColor', mc{mod(ii,size(mc,1))}); % plot the experimental data

    figure(chii_fig)
    scatter(continu_var, medfilt1(analysis.xi), 'Marker', mkr(mod(ii,length(mkr))),...
        'MarkerEdgeColor', mc{mod(ii,size(mc,1))}, 'MarkerFaceColor', mc{mod(ii,size(mc,1))}); % plot the experimental data
%     figure(chii_fig)
%     plot(continu_var(1,1:15:end), analysis.xi(1:15:end), 'LineStyle', '-', 'Marker', mkr(mod(ii,length(mkr))),...
%         'Color', mc{mod(ii,size(mc,1))}); % plot the experimental data

    sim_file = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_4*GHz', temps(ii), theta, phi, gama);
    if Options.RPA == true
        sim_load = strcat("RPA_LiHoF4_", sim_file, ".mat");
        xlabr = 'Re[\chi^{zz}_{RPA}] (meV)';
        xlabi = 'Im[\chi^{zz}_{RPA}] (meV)';
    else
        sim_load = strcat("chi0_LiHoF4_", sim_file, ".mat");
        xlabr = 'Re[\chi^{zz}_{MF}] (meV)';
        xlabi = 'Im[\chi^{zz}_{MF}] (meV)';
    end
    chi_file = dir(fullfile(sim_loc, sim_load)); % use parts unspecified file name
    if isempty(chi_file)
        disp('Cannot match the file, check the name and directory!')
        break
    end
    load([sim_loc chi_file.name],'-mat','freq_total','fields','chi','unit'); % load calculated suscpetibilities
    chiz = squeeze(chi(3,3,:,:));

    % Convert the unit to meV
    if strcmp(unit, 'GHz')
        ConvUnit = 1/Gh2mV;
    elseif strcmp(unit, 'J')
        ConvUnit = meV2J;
    else
        ConvUnit = 1;
    end
    
    [~, fidx] = min(abs(freq_total - fc));
    figure(chir_fig);
    plot(fields, real(chiz(fidx,:))/ConvUnit, 'LineStyle', lin(1), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);

    figure(chii_fig);
    plot(fields, imag(chiz(fidx,:))/ConvUnit, 'LineStyle', lin(1), 'Color', mc{mod(ii,size(mc,1))}, 'lineWidth', 1.5);
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
