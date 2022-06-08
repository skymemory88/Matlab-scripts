%% Make a bar graph of coupling strenght and spin linewidth of electro-nuclear transitions in LiHoF4
clearvars

fc = 3.642; % cavity resonant frequency [GHz]
theta = 0;
phi = 2.3;
n_trans = 6;

sim_temps = [0.11, 0.14, 0.26]; % temperatures
exp_temps = [0.15, 0.18, 0.30]; % temperatures
lgd_sim = string(sim_temps.*1000);
lgd_sim = 'Sim.' + lgd_sim + ' mK';
lgd_exp = string(exp_temps.*1000);
lgd_exp = 'Exp.' + lgd_exp + ' mK';

Options.barfig = 1; % 1: 1D grouped, 2. 2D stacked, 3. 3D detached.
Options.nZee = false; % whether to include nuclear Zeeman term
currentLoc = pwd;
hbar = 1.055E-34; % Reduced Planck constant [J.s]
J2meV = 6.24151e+21; % [meV/J]
Gh2mV = hbar*2*pi*10^9*J2meV; % GHz to meV

m = ['o','s','^','x','<','>'];
mc = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]];
% mc = ['k','m','y','k','g'];
xnames = ["|5\rangle\rightarrow|6\rangle", "|4\rangle\rightarrow|5\rangle", "|3\rangle\rightarrow|4\rangle",...
    "|2\rangle\rightarrow|3\rangle", "|1\rangle\rightarrow|2\rangle", "|0\rangle\rightarrow|1\rangle"];

if Options.nZee == true
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li', 'Ho',...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I'];
else
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li', 'Ho',...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\without Hz_I'];
end

g_sim = zeros(length(exp_temps), n_trans);
for ii = 1:length(sim_temps)
    cd 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4'
    lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', sim_temps(ii), theta, phi, 1),'.mat'];
    eFile = fullfile(Options.location,lname);
    load(eFile,'-mat','eee','fff');
    eee = eee- min(eee,[],2); % Normalize the energy amplitude to the lowest eigen-energy
    eigenE = squeeze(eee) / Gh2mV; % [GHz]
    eigenE = eigenE(:,1:n_trans+1);
    fields = vecnorm(fff); % Magnetic field [T]
    Ediff = diff(eigenE,1,2); % Transition between the nearest neighbouring levels
    sim = S11_simulation('Ho', 'field', sim_temps(ii), 0, 2.3, 9e-5, false);
    for jj = 1:size(Ediff,2)
        [~,idx] = min(abs(Ediff(:,jj)-fc));
        b0(jj) = fields(idx);
        g_sim(ii,jj) = sqrt(sim.Gc2(idx,jj+1,jj));
    end
%     % save the raw data in an ascii file (.txt)
%     data_sim = [reshape(b0,length(b0),1) reshape(g_sim(ii,:),length(g_sim(ii,:)),1)];
%     loc = 'C:\Users\yiyang\Desktop\';
%     save(fullfile(loc,sprintf('Sim_SC239_%umK_gc.txt',sim_temps(ii)*1000)),'data_sim','-ascii');
end

barGrp = figure;
box on
hold on
xticks([1:6])
g_sim = flip(g_sim,2);
if Options.barfig == 1
    bar(g_sim','grouped')
    set(gca,'XTickLabel',xnames)
    xlabel('Transitions')
    ylabel('Coupling strengh (GHz)')
elseif Options.barfig == 2
    [tempData, rIdx] = sort(g_sim,1,'descend');
    for ii = 1:size(g_sim,1)
        figure(barGrp);
        fig = bar(tempData(ii,:), 'stacked','FaceColor','flat');
        for jj = 1:size(tempData(ii,:),2)
            fig.CData(jj,:) = mc(rIdx(ii,jj),:);
        end
    end
    set(gca,'XTickLabel',xnames)
    xlabel('Transitions')
    ylabel('Coupling strengh (GHz)')
elseif Options.barfig == 3
    bar3([1:6], g_sim', 'detached')
    set(gca,'XTickLabel', sim_temps)
    yticks([1:6])
    set(gca,'YTickLabel', xnames)
    xlabel('Temperature (K)')
    ylabel('Transitions')
    zlabel('Coupling strengh (GHz)')
end

B0 = zeros(length(exp_temps), n_trans);
g_exp = zeros(length(exp_temps), n_trans);
g_err = zeros(length(exp_temps), n_trans);
for ii = 1:length(exp_temps)
    cd 'G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Fieldscan'
    [b_temp, g_temp, ge_temp, ~, ~] = batch_plot(exp_temps(ii), false);
    if exp_temps(ii) == 0.15 % For data at 150 mK, retain spots for the 5-6 transition with zeros
        b_temp = [0 b_temp(2:end)];
        g_temp = [0 g_temp(2:end)];
        ge_temp = [0 ge_temp(2:end)];
    else
        b_temp = b_temp(2:end); % Remove the coupling on FM side
        g_temp = g_temp(2:end);
        ge_temp = ge_temp(2:end);
    end

    for jj = 1:length(b_temp)
        B0(ii,jj) = b_temp(jj);
        g_exp(ii,jj) = g_temp(jj);
        g_err(ii,jj) = ge_temp(jj);
    end
end


xx = 1:6;
if Options.barfig == 1
    xx = xx - 0.2;
    for ii = 1:length(exp_temps)
        xp = xx(g_exp(ii,:) ~= 0);
        xx = xx + 0.2;
%         plot(xp, nonzeros(g_exp(ii,:)), 'Marker', m(ii), 'MarkerFaceColor', mc(ii,:), 'MarkerEdgeColor', 'k');
        plot(xp, nonzeros(g_exp(ii,:)), 'Marker', m(ii), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    end
elseif Options.barfig == 2
    for ii = 1:length(exp_temps)
        xp = xx(g_exp(ii,:) ~= 0);
        plot(xp, nonzeros(g_exp(ii,:)), 'Marker', m(ii), 'MarkerFaceColor', mc(ii,:), 'MarkerEdgeColor', 'k');
    end
elseif Options.barfig == 3
    for ii = 1:length(exp_temps)
        xp = xx(g_exp(ii,:) ~= 0);
        [~,xidx] = min(abs(sim_temps-exp_temps(ii)+0.05));
        plot3(xidx*ones(1,length(xp)), xp, nonzeros(g_exp(ii,:)), 'Marker', m(ii),...
            'MarkerFaceColor', mc(ii,:), 'MarkerEdgeColor', 'k', 'linewidth', 2);
        view(140,30);
    end
end

legend([lgd_sim lgd_exp], 'Location', 'northwest')
cd(currentLoc)