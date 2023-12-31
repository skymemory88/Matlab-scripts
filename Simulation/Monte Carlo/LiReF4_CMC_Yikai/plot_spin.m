function plot_spin(init, meas, params, ff, plotOpt, spinOpt)

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
    pOpt = plotOpt(ii);
    if strcmpi(pOpt, 'domain') || strcmpi(pOpt, 'all')
        fig1 = figure;
        sf1 = subplot(1,3,1);
        Sz = squeeze(sign(spin0(:,3,:))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc(:,ff)) );
        color = [cc(:,ff) c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(params.pos(:,1), params.pos(:,2), params.pos(:,3), 35, color, 'filled');
        pbaspect([1 1 1])
        title(sprintf('Initial domain config at B = %.2f T', vecnorm(params.field(:,ff))))
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.04 0.5 0.5]);

        sf2 = subplot(1,3,2);
        Sz = squeeze(sign(spinT(:,3,:))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc(:,ff)) );
        color = [cc(:,ff) c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(params.pos(:,1), params.pos(:,2), params.pos(:,3), 35, color, 'filled');
        pbaspect([1 1 1])
        title(sprintf('Thermalized domain config at B = %.2f T', vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        Sz = squeeze(sign(spin(:,3,:))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc(:,ff)) );
        color = [cc(:,ff) c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(params.pos(:,1), params.pos(:,2), params.pos(:,3), 35, color, 'filled');
        pbaspect([1 1 1])
        title(sprintf('Final domain config at B = %.2f T', vecnorm(params.field(:,ff))))
    end
    if strcmpi(pOpt, 'vector') || strcmpi(pOpt, 'all')
        fig2 = figure;
        sf1 = subplot(1,3,1);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spin0(:,1,1,ff), spin0(:,2,1,ff), spin0(:,3,1,ff), 'k');
        pbaspect([1 1 1])
        title(sprintf('Initial spin config at B = %.2f T', vecnorm(params.field(:,ff))))
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.4 0.5 0.5]);

        sf2 = subplot(1,3,2);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spinT(:,1,1,ff), spinT(:,2,1,ff), spinT(:,3,1,ff), 'b');
        pbaspect([1 1 1])
        title(sprintf('Thermalized spin config at B = %.2f T', vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        quiver3(params.pos(:,1), params.pos(:,2), params.pos(:,3), spin(:,1,1,ff), spin(:,2,1,ff), spin(:,3,1,ff),'r');
        pbaspect([1 1 1])
        title(sprintf('Final spin config at B = %.2f T', vecnorm(params.field(:,ff))))
    end
    if strcmpi(pOpt, 'magnetization') || strcmpi(pOpt, 'all')
        fg3 = figure;
        ax3 = fg3.CurrentAxes;
        plot(vecnorm(params.field,1), squeeze(rms(meas.Mx,1)), 'o');
        hold on
        plot(vecnorm(params.field,1), squeeze(rms(meas.My,1)), 'o');
        plot(vecnorm(params.field,1), squeeze(rms(meas.Mz,1)), 'o');
        legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
        xlabel('Magnetic Field (T)')
        ylabel('$\langle J \rangle$','Interpreter','latex')
        set(ax3, 'FontSize', 12);
    end
end
end