%% Make a bar graph of coupling strenght and spin linewidth of electro-nuclear transitions in LiHoF4
clearvars

Options.figSty = 'Bar1D'; % 1: Bar1D grouped, 2. Bar2D stacked, 3. Bar3D detached, 4. TempFig
Options.nZee = true; % Option to include nuclear Zeeman effect
Options.gamma = true; % Option to plot spin linewidth from fitting
mion = 'Ho'; % Magnetic ion species

fc = 3.642; % cavity resonant frequency [GHz]
Hc = 5; % critical field (T)
theta = 0.0; % magnetic field angle in a-c plane
phi = -13.5; % magnetic field angle in a-b plane
n_trans = 6; % number of excitation spectra to include

sim_temps = [0.1 0.18 0.3]; % simulation temperatures
exp_temps = [0.15 0.18 0.3]; % simulation temperatures
% sim_temps = 0.15; 
% exp_temps = 0.15;
% crossB = [6.502 8.266 9.838 11.43 13.00 14.28]; % manually set field locations for anti-crossings

lgd_sim = string(sim_temps.*1000);
lgd_sim = 'Sim.' + lgd_sim + ' mK';
lgd_exp = string(exp_temps.*1000);
lgd_exp = 'Exp.' + lgd_exp + ' mK';

currentLoc = pwd;
hbar = 1.055E-34; % Reduced Planck constant [J.s]
J2meV = 6.24151e+21; % [meV/J]
Gh2mV = hbar*2*pi*10^9*J2meV; % GHz to meV

mkr = ['o','s','d','^','<','>']; % marker style pool
lin = ['-','--','-.',':']; % line style pool
mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool
% mc = ['k','m','y','k','g']; % marker color pool
xnames = ["|5\rangle\rightarrow|6\rangle", "|4\rangle\rightarrow|5\rangle", "|3\rangle\rightarrow|4\rangle",...
    "|2\rangle\rightarrow|3\rangle", "|1\rangle\rightarrow|2\rangle", "|0\rangle\rightarrow|1\rangle"];

if strcmp(pathsep, ':')
    platform = 'Unix';
    divd = '\';
else
    platform = 'Win';
    divd = '\';
end

% determine the OS envirnment
switch  platform
    case 'Win'
        fpath = ['G:\My Drive\File sharing\PhD program\Research projects\Li', mion, 'F4 project\Data\',...
            'Simulations\Matlab\Susceptibilities\'];
        addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan');
        addpath('');
    case 'Unix'
        fpath = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
            '/PhD program/Research projects/Li', mion, 'F4 project/Data/',...
            'Simulations/Matlab/Susceptibilities/'];
        addpath('/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/File sharing/',...
            'Programming scripts/Matlab/Plot functions/Fieldscan');
        addpath(['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/',...
            'File sharing/Programming scripts/Matlab/Simulation/Mean Field/LiReF4'])
end

% nuclear Zeeman interaction option
if Options.nZee == true
    Options.location = [fpath,'Hz_I=1'];
else
    Options.location = [fpath,'Hz_I=0'];
end

g_sim = zeros(length(exp_temps), n_trans); % coupling strengths
pop = zeros(length(exp_temps), n_trans); % Population factor (P_n - P_n+1)
g2s_sim = zeros(length(exp_temps), 1); % coupling strengths
for ii = 1:length(sim_temps)
    lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', sim_temps(ii), theta, phi, 1),'.mat'];
    eFile = fullfile(Options.location, lname);
    load(eFile,'-mat','eee','fff');
    eee = eee- min(eee,[],2); % Normalize the energy amplitude to the lowest eigen-energy
    eigenE = squeeze(eee) / Gh2mV; % [GHz]
    eigenE = eigenE(:,1:n_trans+1);
    fields = vecnorm(fff); % Magnetic field [T]
    Ediff = diff(eigenE,1,2); % Transition between the nearest neighbouring levels
%     calc = S11_simulation('Ho', 'field', sim_temps(ii), 0, 2.3, 1e-4, false);
%     Gc2 = poermute(calc.Gc2{:},[2,3,1]);
    custom.Mx = true;
    custom.Enorm = true;
    calc = EngyLevels(sim_temps(ii), theta, phi, 7, custom);
    NN = squeeze(calc.pop{:});
    Gc2 = squeeze(calc.Gc2{:});
%     b0 = double.empty(size(Ediff,2),0); % f
%     wDif = double.empty(size(Ediff,2),0); % debug tool: cavity-spin mode crossing checkpoint
    for jj = 1:size(Ediff,2)
%         [~,idx] = min(abs(fields-crossB(jj))); % use manually set field locations
%         [~,idx] = min(abs(Ediff(:,jj)-fc)); % search at the entire field range
%         [~,idx] = min(abs(Ediff(fields <= Hc,jj)-fc)); % search only below critical field (FM state)
        [~,idx] = min(abs(Ediff(fields >= Hc,jj)-fc)); % search only above critical field (PM state)
        idx = length(fields(fields < Hc)) + idx; % add back the truncated length
%         b0(jj) = fields(idx); % for debugging
        g_sim(ii,jj) = sqrt(Gc2(jj+1, jj, idx));
        pop(ii,jj) = squeeze(NN(jj+1, jj, idx)); % find the appropriate excitation mode
    end
    g2s_sim(ii) = sqrt( sum(squeeze(g_sim(ii,:)).^2,2) );
%     % save the raw data in an ascii file (.txt)
%     data_sim = [reshape(b0,length(b0),1) reshape(g_sim(ii,:),length(g_sim(ii,:)),1)];
%     loc = 'C:\Users\yiyang\Desktop\';
%     save(fullfile(loc,sprintf('Sim_SC239_%umK_gc.txt',sim_temps(ii)*1000)),'data_sim','-ascii');
end

B0 = zeros(length(exp_temps), n_trans);
g_exp = zeros(length(exp_temps), n_trans);
g_err = zeros(length(exp_temps), n_trans);
gma_exp = zeros(length(exp_temps), n_trans);
gFM_exp = zeros(length(exp_temps), 0);
gmaFM_exp = zeros(length(exp_temps), 0);
for ii = 1:length(exp_temps)
    [b_temp, g_temp, ge_temp, gma_temp, ~] = batch_plot(exp_temps(ii), false, 'on');
    gFM_exp(ii) = g_temp(1); % record the coupling strength on the FM side
    gmaFM_exp(ii) = gma_temp(1);

    b_temp = b_temp(2:end); % Remove the coupling on FM side
    g_temp = g_temp(2:end);
    ge_temp = ge_temp(2:end);
    gma_temp = gma_temp(2:end);

    for jj = 1:length(b_temp)
        B0(ii,jj) = b_temp(jj);
        g_exp(ii,jj) = g_temp(jj);
        g_err(ii,jj) = ge_temp(jj);
        gma_exp(ii,jj) = gma_temp(jj);
    end
end


Grp = figure;
box on
hold on
xticks([1:6])
g_sim = flip(g_sim,2);
pop = flip(pop,2);
xx = 1:6;
switch Options.figSty
    case 'Bar1D'
        bar(g_sim','grouped')
        xx = xx - 0.05*length(exp_temps);
        for ii = 1:length(exp_temps)
            xp = xx(g_exp(ii,:) ~= 0);
            xx = xx + 0.15;
            plot(xp, nonzeros(g_exp(ii,:)), 'Marker', mkr(ii), 'MarkerFaceColor', mc{ii}, 'MarkerSize', 10,...
                'MarkerEdgeColor', 'k', 'Color', 'k');
        end
        set(gca,'XTickLabel',xnames)
        xlabel('Transitions')
        ylabel('Coupling strengh (GHz)')
        legend([lgd_sim lgd_exp], 'Location', 'northwest')
    case 'Bar2D'
        [tempData, rIdx] = sort(g_sim,1,'descend');
        for ii = 1:size(g_sim,1)
            figure(Grp);
            fig = bar(tempData(ii,:), 'stacked','FaceColor','flat');
            for jj = 1:size(tempData(ii,:),2)
                fig.CData(jj,:) = mc{rIdx(ii,jj)};
            end
        end
        for ii = 1:length(exp_temps)
            xp = xx(g_exp(ii,:) ~= 0);
            plot(xp, nonzeros(g_exp(ii,:)), 'Marker', mkr(ii), 'MarkerFaceColor', mc{ii}, 'MarkerEdgeColor', 'k');
        end
        set(gca,'XTickLabel',xnames)
        xlabel('Transitions')
        ylabel('Coupling strengh (GHz)')
        legend([lgd_sim lgd_exp], 'Location', 'northwest')
    case 'Bar3D'
        bar3([1:6], g_sim', 'detached')
        for ii = 1:length(exp_temps)
            xp = xx(g_exp(ii,:) ~= 0);
            [~,xidx] = min(abs(sim_temps-exp_temps(ii)+0.05));
            plot3(xidx*ones(1,length(xp)), xp, nonzeros(g_exp(ii,:)), 'Marker', mkr(ii),...
                'MarkerFaceColor', mc{ii}, 'MarkerEdgeColor', 'k', 'linewidth', 2);
            view(140,30);
        end
        set(gca,'XTickLabel', sim_temps)
        yticks([1:6])
        set(gca,'YTickLabel', xnames)
        xlabel('Temperature (K)')
        ylabel('Transitions')
        zlabel('Coupling strengh (GHz)')
        legend([lgd_sim lgd_exp], 'Location', 'northwest')
    case 'TempFig'
        plot(sim_temps, g_sim, '-');
        plot(sim_temps, g2s_sim, '--k');        
        for ii = 1:length(exp_temps)
            fcol = mc{g_exp(ii,:) ~= 0};
            scatter(exp_temps(ii), nonzeros(g_exp(ii,:)), 30, fcol, 'filled','MarkerEdgeColor', 'k');
%             extsn = plot([0.053 exp_temps(ii)], [nonzeros(g_exp(ii,:)) nonzeros(g_exp(ii,:))]', 'lineStyle', lin(ii));
%             for jj = 1:size(extsn)
%                 set(extsn(jj),'color',fcol(jj,:))
%             end
        end
        xlabel('Temperature (K)')
        ylabel('Coupling strengh (GHz)')
end

if Options.gamma == true
    xx = 1:6;
    gfig = figure;
    hold on
    box on
    cfig = figure;
    hold on
    box on
    for ii = 1:length(exp_temps)
        xp = xx(gma_exp(ii,:) ~= 0);
        figure(gfig)
        plot(xp, nonzeros(gma_exp(ii,:)), 'Marker', mkr(ii), 'MarkerFaceColor', mc{ii}, 'MarkerSize', 6,...
            'MarkerEdgeColor', 'k');
        figure(cfig)
        plot(xp, nonzeros(g_exp(ii,:)).^2./nonzeros(gma_exp(ii,:))/5e-4, 'Marker', mkr(ii), 'MarkerFaceColor', mc{ii}, 'MarkerSize', 6,...
            'MarkerEdgeColor', 'k');
    end
    figure(gfig)
    ylabel('Spin line width (GHz)')
    set(gca,'XTickLabel',xnames)
    legend([lgd_exp], 'Location', 'northwest')
    figure(cfig)
    ylabel('Cooperativity')
    set(gca,'XTickLabel',xnames)
    legend([lgd_exp], 'Location', 'northwest')
end
cd(currentLoc)