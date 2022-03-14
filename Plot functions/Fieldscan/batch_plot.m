function [B0, gcs, g_ratio] = batch_plot(temp)

Options.cmap = false;
Options.data = false;
Options.ebar = false;
Options.mean = true;
Options.conden = false;

if temp == 0.15
    aux_files(1) = "2021_07_0003_0.148K_1.01-4.89T";
%     aux_files(2) = "2021_07_0003_0.142K_7.51-8.49T";
    aux_files(2) = "2021_07_0003_0.148K_6.80-8.79T";
    aux_files(3) = "2021_07_0003_0.142K_8.80-10.19T";
    aux_files(4) = "2021_07_0003_0.142K_10.31-12.10T";
    aux_files(5) = "2021_07_0003_0.142K_12.11-13.49T";
%     aux_files(6) = "2021_07_0003_0.142K_12.80-17.00T";
    aux_files(6) = "2021_07_0003_0.142K_13.01-17.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.04\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_150mK_-25dBm_0dB.fig';
elseif temp == 0.18
    aux_files(1) = "2021_07_0004_0.167K_2.00-4.50T";
    aux_files(2) = "2021_07_0004_0.167K_5.51-7.00T";
    aux_files(3) = "2021_07_0004_0.167K_7.21-9.00T";
    aux_files(4) = "2021_07_0004_0.167K_8.51-10.70T";
    aux_files(5) = "2021_07_0004_0.167K_10.41-12.59T";
    aux_files(6) = "2021_07_0004_0.167K_12.11-13.79T";
%     aux_files(6) = "2021_07_0004_0.167K_12.00-13.79T";
    aux_files(7) = "2021_07_0004_0.167K_12.50-15.99T";
%     aux_files(7) = "2021_07_0004_0.167K_12.71-17.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.05\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_180mK_-25dBm_0dB.fig';
end

B0 = double.empty(length(aux_files),0);
gcs = double.empty(length(aux_files),0);
for ii = 1:length(aux_files)
    load_Obj = fullfile(location,strcat(aux_files(ii),".mat"));
    load(load_Obj,'H0','w0','w0_ci','gc','gc_ci','gma','gma_ci');
    B0(ii) = mean(H0);
    if ii == 1
        gplot = figure;
        hold on
        xlabel('Magnetic Field (T)')
        ylabel('Coupling strength (GHz)')
    else
        figure(gplot);
        if Options.conden == true; yyaxis left; end
    end
    if Options.data == true
        if Options.ebar == true
            errorbar(H0,gc,abs(gc_ci(:,1)-gc_ci(:,2)),'ok');
        else
            plot(H0,gc,'.k');
        end
    end
    if Options.mean == true
        plot(mean(H0),mean(gc),'sb');
        gcs(ii) = mean(gc);
    end

    if ii == 1
        if Options.conden == true
            figure(gplot);
            yyaxis right
            hold on
            ylabel('Spin linewidth (GHz)')
        else
            rplot = figure;
            hold on
            xlabel('Magnetic Field (T)')
            ylabel('Spin linewidth (GHz)')
        end
    else
        if Options.conden == true
            figure(gplot);
            yyaxis right
        else
            figure(rplot);
        end
    end
    if Options.data == true
        if Options.ebar == true
            errorbar(H0,gma,abs(gma_ci(:,1)-gma_ci(:,2)),'ok');
        else
            plot(H0,gma,'.k');
        end
    end
    if Options.mean == true; plot(mean(H0),mean(gma),'sr'); end

    if Options.cmap == true
        if ii == 1
            cmap = open(strcat(location,cmapfig));
            figure(cmap)
            colormap('bone')
            hold on
            plot(H0,w0,'.r');
        else
            figure(cmap);
            plot(H0,w0,'.r');
        end
%         if Options.ebar == true; plot(H0,w0_ci(:,1),'-r',H0,w0_ci(:,2),'-r'); end
    end
end
if Options.mean == true
    rplot = figure;
    g_ratio = zeros(length(aux_files)-1,1);
    for ii = 2:length(aux_files)-1 % excluding the anti-crossing on the FM side
    figure(rplot)
    hold on
    g_ratio(ii) = gcs(ii)./gcs(ii+1);
    plot(ii, g_ratio(ii),'o')
    xlabel('Eigen-states')
    ylabel('Gc_i/Gc_j ratio (arb.)')
    end
end

% plot theoretical values of g-ratios
wc = 3.645;
hbar = 1.055E-34; % Reduced Planck constant [J.s]
meV2J = 1.602217e-22; % [J/meV]
kB = 0.08617; % [meV/K]
hhbar = hbar/meV2J;
hhbar = hhbar*1e9;
tt = 0:0.001:0.5;
figure;
plot(tt,exp(-hhbar*wc/kB./tt),'.k');
xlabel('Temperature (K)')
ylabel('exp(-\betah\omega_c)')
hold on
% plot(temp,exp(-hhbar*3.645/kB./temp),'or');
legend('exp(-\beta h\omega_c)','150 mK','180 mK')

end