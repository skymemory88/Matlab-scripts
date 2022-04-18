function [B0, gcs, gc_err, gama, gama_err] = batch_plot(temp, fig_opt)
% batch retrieve fitted parameters from multiple measurements

Options.figure = fig_opt;
Options.colmap = false;
Options.data = true;
Options.ebar = false;
Options.mean = true;
Options.conden = false;

[aux_files, location, cmapfig] = fit_list(temp);

B0 = double.empty(length(aux_files),0);
gcs = double.empty(length(aux_files),0);
gc_err = double.empty(length(aux_files),0);
gama = double.empty(length(aux_files),0);
gama_err = double.empty(length(aux_files),0);
for ii = 1:length(aux_files)
    load_Obj = fullfile(location,strcat(aux_files(ii),".mat"));
    load(load_Obj,'H0','w0','w0_ci','gc','gc_ci','gma','gma_ci');
    B0(ii) = mean(H0);
    gc_err(ii) = mae(abs(gc_ci(:,1)-gc_ci(:,2)));
    gcs(ii) = mean(gc);
    gama(ii) = mean(gma);
    gama_err(ii) = mae(abs(gma_ci(:,1)-gma_ci(:,2)));
%     % save the raw data in an ascii file (.txt)
%     g_fit = [reshape(H0,length(H0),1) reshape(gc,length(gc),1) reshape(gma,length(gma),1)];
%     loc = 'C:\Users\yiyang\Desktop\';
%     save(fullfile(loc,sprintf('SC239_%umK_%u_gFit.txt',temp*1000,ii)),'g_fit','-ascii');

    if Options.figure == true
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
                errorbar(H0,gc,abs(gc_ci(:,1)-gc_ci(:,2)),'.k');
            else
                plot(H0,gc,'.k');
            end
        end
        if Options.mean == true
            if Options.ebar == true
    %             errorbar(mean(H0), mean(gc), gc_err(ii),'or');
                errorbar(mean(H0), mean(gc), std(gc)/sqrt(length(gc)),'or');
    %             errorbar(mean(H0), mean(gc), std(gc),'or');
            else
                plot(mean(H0),mean(gc),'or');
            end
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
                errorbar(H0, gma, abs(gma_ci(:,1)-gma_ci(:,2)),'.k');
            else
                plot(H0, gma, '.k');
            end
        end
        if Options.mean == true
            if Options.ebar == true
%                 errorbar(mean(H0), mean(gma), gama_err(ii),'or');
                errorbar(mean(H0), mean(gma), std(gma)/sqrt(length(gma)),'or');
%                 errorbar(mean(H0), mean(gma), std(gma),'or');
            else
                plot(mean(H0),mean(gma),'or');
            end
        end

        if Options.colmap == true
            if ii == 1
                colmap = open(strcat(location,cmapfig));
                figure(colmap)
                hold on
                plot(H0,w0,'.r');
            else
                figure(colmap);
                plot(H0,w0,'.r');
            end
%             cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
%                 [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
%                 [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
%             cmap = flip(cmap,1);
%             colormap(cmap)
%             colormap('bone')
            colormap('jet')
%             if Options.ebar == true; plot(H0,w0_ci(:,1),'-r',H0,w0_ci(:,2),'-r'); end
        end
    end
end
if Options.figure == true && Options.mean == true
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

% % plot theoretical values of g-ratios
% wc = 3.645;
% hbar = 1.055E-34; % Reduced Planck constant [J.s]
% meV2J = 1.602217e-22; % [J/meV]
% kB = 0.08617; % [meV/K]
% hhbar = hbar/meV2J;
% hhbar = hhbar*1e9;
% tt = 0:0.001:0.5;
% figure;
% plot(tt,exp(-hhbar*wc/kB./tt),'.k');
% xlabel('Temperature (K)')
% ylabel('exp(-\betah\omega_c)')
% hold on
% % plot(temp,exp(-hhbar*3.645/kB./temp),'or');
% legend('exp(-\beta h\omega_c)','150 mK','180 mK')

end