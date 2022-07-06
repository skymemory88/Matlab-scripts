Options.RPA  = true;
Options.nZee = false; % Nuclear Zeeman interaction
Options.sample = 239; % Sample to use (SC###)

temp = [0.1 0.3 0.8 1]; % temperatures
% temp = 1;
theta = 0.0; % a-c plane rotation
phi = 2.3; % a-b plane rotation
freq_l = 4.70; % frequency range lower limit [GHz]
freq_h = 4.745; % frequency range upper limit [GHz]
gama = 1e-4; % default spin linewidth

hbar = 1.05457E-34; % Reduced Planck constant [J.s/rad]
meV2J = 1.60218e-22; % [J/meV]
Gh2mV = hbar * 2*pi / meV2J * 10^9; % [meV/GHz]

% set experimental data loading parameters
exp_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC',...
    num2str(Options.sample),'\'];

dat_lst = readtable(fullfile(exp_loc, "off_res_list.xlsx"));
exp_fold = table2array(dat_lst(:,1));
exp_files = table2array(dat_lst(:,2));
exp_temps = table2array(dat_lst(:,3));

% set calculated data loading parameters
if Options.nZee == true
    sim_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\',...
        'Susceptibilities\with Hz_I\backup\'];
else
    sim_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab\',...
        'Susceptibilities\without Hz_I\backup\'];
end

% plot the data
chi_fig = figure;
box on
hold on
for ii = 1:length(temp)
    [~, idx] = min(abs(exp_temps-temp(ii))); % find the closest experimental temperature
    exp_file = fullfile([exp_loc, exp_fold{idx},'\'], exp_files{idx}); % locate the experimental data
    load([exp_file, '_interp'], '-mat', 'analysis', 'continu_var');
    fc = nonzeros(analysis.wc);
    fc = fc(1); % pick zero-field resonant frequency as the reference frequency
    figure(chi_fig)
    plot(continu_var(1,:), analysis.xr, '.r'); % plot the experimental data

    sim_file = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_4*GHz', temp(ii), theta, phi, gama);
    if Options.RPA == true
        sim_load = strcat('RPA_LiHoF4_', sim_file, '.mat');
        xlab = 'Re[\chi^{zz}_{RPA}] (meV)';
    else
        sim_load = strcat('chi0_LiHoF4_', sim_file, '.mat');
        xlab = 'Re[\chi^{zz}_{MF}] (meV)';
    end
    chi_file = dir(fullfile(sim_loc, sim_load)); % use parts unspecified file name
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
    figure(chi_fig)
    plot(fields, real(chiz(fidx,:))/ConvUnit,'-k','lineWidth',1.5);
end
xlim([0 9])
xlabel('Magnetick Field (T)');
ylabel(xlab);
