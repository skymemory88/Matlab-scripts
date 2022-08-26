%% isolated frequency sweep plotter
clearvars
B0 = [0.1 4.23 7.0]; % field location of frequency cuts to plot/fit
bidx = uint32.empty(length(B0),0);
% load the processed data
Options.plot = 'compressed'; % plot the frequency sweeps onto the same or individual plots
Options.yScale = 'dB'; % 1. 'lin': linear. 2. 'dB': log scale
Options.fit = true; % Option to fit the data
    dType = 'off'; % specify the data type (on/off resonance) for fitting function selection
Options.load_method = 'manual'; % 1: auto, 2: manual, 3: skip (if already exist in workspace). 
    temp = 0.500; % temperature (K) for locating the file in 'auto' mode
    sample = 239; % sample SC number for locating the file in 'auto' mode

mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool

switch Options.load_method
    case 'auto'
        exp_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC',...
            num2str(sample), '\'];
        dat_lst = readtable(fullfile(exp_loc, sprintf("SC%u_off_list.xlsx", sample)));
        folder = table2array(dat_lst(:,1));
        efiles = table2array(dat_lst(:,2));
        etemps = table2array(dat_lst(:,3));
        [~, idx] = min(abs(etemps-temp)); % find the closest experimental temperature
        exp_file = fullfile([exp_loc, folder{idx},'\'], efiles{idx}); % locate the experimental data
        load([exp_file, '_interp'], '-mat', 'S11', 'freq', 'continu_var', 'analysis');
    case 'manual'
        fileloc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\',...
            'File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
            'SC107 (4x5x2mm)\19.05.2019'];
%         fileloc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\',...
%             'SC239\2021.07.16'];
        loadname = '2019_05_0026';
        LoadObj = fullfile(fileloc,loadname);
        load([LoadObj, '_interp'], '-mat', 'S11', 'freq', 'continu_var', 'analysis');
        for ii = 1:size(S11,2)
            S11(:,ii) = S11(:,ii)./max(S11(:,ii));
%             S11(:,ii) = S11(:,ii)./0.9;
        end
end

% define colormap
cmap = unique([[0 0 0];[zeros(200,1),linspace(0,0.5,200)',linspace(0.4,1,200)'];...
    [linspace(0,1,200)',0.5*ones(200,1),linspace(1,0,200)'];...
    [ones(200,1),linspace(0.5,1,200)',linspace(0,1,200)']],'rows');
cmap = flip(cmap,1);

sfig = figure('Colormap', cmap, 'Position',[200 200 800 650]);
box on
hold on
cplot = pcolor(continu_var, freq, mag2db(S11));
cplot.EdgeColor = 'none';
c = colorbar('Location','northoutside');
pos = get(gca,'Position');
set(c,'Position',[pos(1) pos(2)+pos(4) pos(3)-0.1 0.05]);
set(gca,'Position',[pos(1) pos(2) pos(3) pos(4)-0.01])
c.Label.String = '|S11|';
c.Label.Position = [pos(3)+1.5 0.2];
set(gca, 'FontSize', 14)
xlabel('Magnetic Field (T)')
ylabel('Frequency (GHz)')

if Options.fit == true
    [~, idx0] = min(continu_var(1,:));
    weight(:,:,1) = abs(1./S11(:,:));      % Weight function option 1
    % weight = 1./gradientweight(S11);       % Weight function option 2
    % weight = ones(size(S11));                % Weight function option 3 (uniform)
    weight(isinf(weight)) = 1e4;           % Remove infinities in weight function
    param  =  [1e-3    1e-4    0     0    freq(S11(:,idx0)==min(S11(:,idx0)),1)]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
    bound_l = [ 0       0      0     0    min(freq(:,idx0))]; % lower bound of fitting parameters
    bound_h = [Inf     Inf     0     0    max(freq(:,idx0))]; % upper bound of fitting parameters
    fit0 = iptoptx(freq(:,idx0), S11(:,idx0), continu_var(1,idx0), param, bound_l, bound_h, weight(:,idx0), false);
    kpe0 = fit0.kpe;
    kpi0 = fit0.kpi;
    f_init = fit0.wc;
end

S11_fit = cell(length(B0),1);
for ii = 1:length(B0)
    [~, bidx(ii)] = min(abs(B0(ii)-continu_var(1,:)));
    B0(ii) = continu_var(1,bidx(ii));
    figure(sfig)
    plot([B0(ii) B0(ii)], [min(freq(:,ii)) max(freq(:,ii))], '--w', 'LineWidth', 2)

    if Options.fit == true
        % Fit using input-output formalism
        param  =  [kpe0  kpi0     0      0    f_init]; % param = {'kpe', 'kpi', 'Re[chi]', 'Im[chi]', 'wc'}
        bound_l = [kpe0-0.1  kpi0   -Inf   -Inf   f_init]; % lower bound of fitting parameters
        bound_h = [kpe0+0.1  kpi0    Inf    Inf   f_init]; % upper bound of fitting parameters
        S11_fit{ii} = iptoptx(freq(:,bidx(ii)), S11(:,bidx(ii)), B0(ii), param, bound_l, bound_h, weight(:,bidx(ii)), false);
    end
end

switch Options.plot
    case 'compressed'
        fig = figure('Position',[1000 200 350 650]);
        box on
        hold on
        for ii = 1:length(B0)
            switch Options.yScale
                case 'dB'
                    yDat = mag2db(S11(:,bidx(ii)));
                    yFit = mag2db(feval(S11_fit{ii},freq(:,bidx(ii))));
                case 'lin'
                    yDat = S11(:,bidx(ii));
                    yFit = feval(S11_fit{ii},freq(:,bidx(ii)));
            end
            scatter( freq(:,bidx(ii)), yDat, 'o', 'filled');
%             plot(freq(:,bidx(ii)), yFit, '-', 'Color', mc{2});
            plot(freq(:,bidx(ii)), yFit, '--k', 'LineWidth', 2);
        end
        set(gca, 'FontSize', 14)
        xlabel('Frequency (GHz)')
        ylabel('|S11| (dB)')
        camroll(90)
        set(gca,'FontSize',14)
    case 'individual'
        fig = gobjects(length(B0),2);
        for ii = 1:length(B0)
            fig(ii) = figure('Position',[1000 200 350 650]);
            box on
            hold on
            switch Options.yScale
                case 'dB'
                    yDat = mag2db(S11(:,bidx(ii)));
                    yFit = mag2db(feval(S11_fit{ii},freq(:,bidx(ii))));
                case 'lin'
                    yDat = S11(:,bidx(ii));
                    yFit = feval(S11_fit{ii},freq(:,bidx(ii)));
            end
            scatter( freq(:,bidx(ii)), yDat, 'o', 'filled');
            plot(freq(:,bidx(ii)), yFit, '-', 'Color', mc{2});
            set(gca, 'FontSize', 14)
            xlabel('Frequency (GHz)')
            ylabel('|S11| (dB)')
            camroll(90)
            set(gca,'FontSize',14)
        end
end

figure;
plot(continu_var(1,:),analysis.xr,'.');
hold on
ylabel('Re[\chi]')
xlabel('Magnetic Field (T)');
yyaxis right
plot(continu_var(1,:),analysis.xi,'.');
ylabel('Im[\chi]')
set(gca,'FontSize',14)