%% isolated frequency sweep plotter

B0 = [0.1 4.916 5.89]; % field location of frequency cuts to plot/fit
bidx = uint32.empty(length(B0),0);
% load the processed data
Option.plot = 'compressed'; % plot the frequency cuts onto the same or individual plots
Option.load_method = 'manual'; % 1: auto, 2: manual, 3: skip (if already exist in workspace). 
    temp = 0.130; % temperature (K) for locating the file in 'auto' mode
    sample = 239; % sample SC number for locating the file in 'auto' mode

switch Option.load_method
    case 'auto'
        exp_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC',...
            num2str(Options.sample), '\'];
        dat_lst = readtable(fullfile(exp_loc, sprintf("SC%u_off_list.xlsx", sample)));
        folder = table2array(dat_lst(:,1));
        efiles = table2array(dat_lst(:,2));
        etemps = table2array(dat_lst(:,3));
        [~, idx] = min(abs(etemps-temp)); % find the closest experimental temperature
        exp_file = fullfile([exp_loc, folder{idx},'\'], efiles{idx}); % locate the experimental data
        load([exp_file, '_interp'], '-mat', 'S11', 'freq', 'continu_var', 'analysis');
    case 'manual'
        fileloc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
            'SC239\2021.07.16'];
        loadname = '2021_07_0015';
        LoadObj = fullfile(fileloc,loadname);
        load([LoadObj, '_interp'], '-mat', 'S11', 'freq', 'continu_var', 'analysis');
end

% define colormap
cmap = unique([[0 0 0];[zeros(200,1),linspace(0,0.5,200)',linspace(0.4,1,200)'];...
    [linspace(0,1,200)',0.5*ones(200,1),linspace(1,0,200)'];...
    [ones(200,1),linspace(0.5,1,200)',linspace(0,1,200)']],'rows');
cmap = flip(cmap,1);

sfig = figure('Colormap', cmap);
box on
hold on
cplot = pcolor(continu_var, freq, mag2db(S11));
cplot.EdgeColor = 'none';
xlabel('Magnetic Field (T)')
ylabel('Frequency (GHz)')

for ii = 1:length(B0)
    [~, bidx(ii)] = min(abs(B0(ii)-continu_var(1,:)));
    B0(ii) = continu_var(1,bidx(ii));
    figure(sfig)
    plot([B0(ii) B0(ii)], [min(freq(:,ii)) max(freq(:,ii))], '--w')
end

switch Option.plot
    case 'compressed'
        fig = figure('Position',[300 50 350 650]);
        box on
        hold on
        for ii = 1:length(B0)
            scatter( freq(:,bidx(ii)), mag2db(S11(:,bidx(ii))), 'o', 'filled');
        end
        xlabel('Frequency (GHz)')
        ylabel('|S11| (dB)')
        camroll(90)
        set(gca,'FontSize',14)
    case 'individual'
        fig = gobjects(length(B0),2);
        for ii = 1:length(B0)
            fig(ii) = figure('Position',[300 50 350 650]);
            box on
            hold on
            scatter( freq(:,bidx(ii)), mag2db(S11(:,bidx(ii))), 'o', 'filled');
            xlabel('Frequency (GHz)')
            ylabel('|S11| (dB)')
            camroll(90)
            set(gca,'FontSize',14)
        end
end