function plot_spin(init, meas, params, ff, plotOpt, spinOpt)
%%
% plotOpt: 'domain', 'thermalization', 'vector', 'B_mag', 'T_mag', 'all'
% spinOpt: 'electron', 'nuclear'
%%
FntSz = 12; % Font Size
colorOption = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]};

if strcmpi(spinOpt, 'electron')
    spin0 = init.eSpin0;
    spinT = init.eSpinT;
    spin = meas.eSpin;
elseif strcmpi(spinOpt, 'nuclear')
    spin0 = init.nSpin0;
    spinT = init.nSpinT;
    spin = meas.nSpin;    
end

for ii = 1:length(plotOpt)
    pOpt = plotOpt{ii};
    if isstring(pOpt); pOpt = convertStringsToChars(pOpt); end
    switch pOpt(1)
        case 'x'
            jj = 1;
        case 'y'
            jj = 2;
        case 'z'
            jj = 3;
    end
    if strcmpi(pOpt(2:end), 'domain') || strcmpi(pOpt, 'all')
        fig1 = figure;
        sf1 = subplot(1,3,1);
        filt = find(spin0(:,jj,ff)); % take non-zero elements
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spin0(filt,jj,ff))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Initial spin-', pOpt(1), ' domain config at B = %.2f T'], vecnorm(params.field(:,ff))))
        set(fig1, 'Units', 'Normalized', 'OuterPosition', [0.5 0.04 0.5 0.5]);

        sf2 = subplot(1,3,2);
        filt = find(spin0(:,jj,ff));
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spinT(filt,jj,ff))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Thermalized spin-', pOpt(1), ' domain config at B = %.2f T'], vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        filt = find(spin0(:,jj,ff));
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spin(filt,jj,ff))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Final spin-', pOpt(1), ' domain config at B = %.2f T'], vecnorm(params.field(:,ff))))
    end
    if strcmpi(pOpt, 'vector') || strcmpi(pOpt, 'all')
        fig2 = figure;
        sf1 = subplot(1,3,1);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spin0(:,1,1,ff), spin0(:,2,1,ff), spin0(:,3,1,ff), 'k');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Initial spin config at B = %.2f T', vecnorm(params.field(:,ff))))
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.4 0.5 0.5]);

        sf2 = subplot(1,3,2);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spinT(:,1,1,ff), spinT(:,2,1,ff), spinT(:,3,1,ff), 'b');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Thermalized spin config at B = %.2f T', vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spin(:,1,1,ff), spin(:,2,1,ff), spin(:,3,1,ff),'r');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Final spin config at B = %.2f T', vecnorm(params.field(:,ff))))
    end
    if strcmpi(pOpt, 'B_mag') || strcmpi(pOpt, 'all')
        fg3 = figure;
        ax3 = fg3.CurrentAxes;
        hold on
        plot(vecnorm(params.field,1), squeeze(mean(meas.Mx,1)), 'o', 'Color', colorOption{1});
        hold on
        plot(vecnorm(params.field,1), squeeze(mean(meas.My,1)), 'o', 'Color', colorOption{2});
        plot(vecnorm(params.field,1), squeeze(mean(meas.Mz,1)), 'o', 'Color', colorOption{3});
        legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
        xlabel('Magnetic Field (T)')
        ylabel('$\langle J \rangle$','Interpreter','latex')
        set(ax3, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'T_mag') || strcmpi(pOpt, 'all')
        fg4 = figure;
        ax4 = fg4.CurrentAxes;
        hold on
        plot(params.temp, squeeze(rms(meas.Mx,1)), 'o', 'Color', colorOption{1});
        hold on
        plot(params.temp, squeeze(rms(meas.My,1)), 'o', 'Color', colorOption{2});
        plot(params.temp, squeeze(rms(meas.Mz,1)), 'o', 'Color', colorOption{3});
        legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
        xlabel('Temperature (K)')
        ylabel('$\langle J \rangle$','Interpreter','latex')
        set(ax4, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'Cv_T') || strcmpi(pOpt, 'all')
        fg5 = figure;
        ax5 = fg5.CurrentAxes;
        E2m = sum(squeeze(meas.E2m), 1); % <E^2>
        Em2 = sum(squeeze(meas.Em).^2, 1); % <E>^2
        hold on
        plot(params.temp, (E2m - Em2)./params.temp, 'o', 'Color', colorOption{1});
        xlabel('Magnetic Field (T)')
        ylabel('$C_v(T)$','Interpreter','latex')
        set(ax5, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'Cv_M') || strcmpi(pOpt, 'all')
        fg6 = figure;
        ax6 = fg6.CurrentAxes;
        E2m = sum(squeeze(meas.E2m), 1); % <E^2>
        Em2 = sum(squeeze(meas.Em).^2, 1); % <E>^2
        hold on
        plot(vecnorm(params.field,1), (E2m - Em2)./vecnorm(params.field,1), 'o', 'Color', colorOption{1});
        xlabel('Magnetic Field (T)')
        ylabel('$C_v(B)$','Interpreter','latex')
        set(ax6, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'thermalization') || strcmpi(pOpt, 'all')
        fg4 = figure;
        ax4 = fg4.CurrentAxes;
        lgd = [];
        hold on;
        for tt = 1:length(params.temp)
            for ff = 1:size(params.field, 2)
                % Plot with specific color/line style and store the handle
                dEt = squeeze(init.dEt(:, tt, ff));
                cutoff = find(dEt,1,'last');
                pl(tt, ff) = plot(1:cutoff, dEt(1:cutoff), 'LineWidth', 1.5);
                lgd{end+1} = ['Temp: ' num2str(params.temp(tt)) ', Field: ' num2str(ff)];
            end
        end
        xlabel('Iteration', 'Interpreter', 'latex');
        ylabel('$\Delta E$', 'Interpreter', 'latex');
        title('Energy Change Evolution over Iterations', 'Interpreter', 'latex');
        set(ax4, 'FontSize', 12);
        grid on;
        legend(lgd, 'Location', 'best');

        fg5 = figure;
        ax5 = fg5.CurrentAxes;
        lgd1 = [];
        hold on;
        for tt = 1:length(params.temp)
            for ff = 1:size(params.field, 2)
                % Plot with specific color/line style and store the handle
                Etot = squeeze(init.Etot(:, tt, ff));
                cutoff = find(Etot,1,'last');
                pl(tt, ff) = plot(1:cutoff, Etot(1:cutoff), 'LineWidth', 1.5);
                lgd1{end+1} = ['Temp: ' num2str(params.temp(tt)) ', Field: ' num2str(ff)];
            end
        end
        xlabel('Iteration', 'Interpreter', 'latex');
        ylabel('$E_{tot}$', 'Interpreter', 'latex');
        title('Total Evolution Over Iterations', 'Interpreter', 'latex');
        set(ax5, 'FontSize', 12);
        grid on;
        legend(lgd1, 'Location', 'best');
    end
end
end