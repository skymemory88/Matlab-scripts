function plot_spin(ion, init, meas, params, TempField, plotOpt, spinOpt)
%%
% plotOpt: 'domain', 'thermalization', 'vector', 'B_mag', 'T_mag', 'all'
% spinOpt: 'electron', 'nuclear'
%%
const.kB = 1.3806e-23; % [J/K]
const.J2meV = 6.24151e+21; % [mev/J]
kB = const.kB * const.J2meV; % [meV]
FntSz = 12; % Font Size
colorOption = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]};

temp = TempField(1);
field = TempField(2);
[~, tt] = min(abs(params.temp - temp));
[~, ff] = min(abs(vecnorm(params.field, 2, 1) - field));
if length(params.temp) == 1 && size(params.field, 2) > 1
    idx = ff;
elseif  length(params.temp) > 1 && size(params.field, 2) == 1
    idx = tt;
else
    idx = [tt, ff];
end

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
        filt = find(spin0(:,jj,idx)); % take non-zero elements
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spin0(filt,jj,idx))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Initial spin-', pOpt(1),...
            ' domain config at T = %.2f K, B = %.2f T'], params.temp(tt), vecnorm(params.field(:,ff))))

        sf2 = subplot(1,3,2);
        filt = find(spin0(:,jj,idx));
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spinT(filt,jj,idx))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Thermalized spin-', pOpt(1),...
            ' domain config at T = %.2f K, B = %.2f T'], params.temp(tt), vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        filt = find(spin0(:,jj,idx));
        pos = params.pos(filt, :);
        Sz = squeeze(sign(spin(filt,jj,idx))); % take signs of the spin configuration
        cc = Sz+1; % +1 --> +2; -1 --> 0
        cc = cc./2; % +2 --> +1; 0 --> 0
        c2 = zeros( size(cc) );
        color = [cc c2 c2]; % red: [1 0 0], black: [0 0 0]
        scatter3(pos(:,1), pos(:,2), pos(:,3), 35, color, 'filled');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf(['Final spin-', pOpt(1),...
            ' domain config at T = %.2f K, B = %.2f T'], params.temp(tt), vecnorm(params.field(:,ff))))
        
        set(fig1, 'Units', 'Normalized', 'OuterPosition', [0.5 0.04 0.5 0.5]);
    end
    if strcmpi(pOpt, 'vector') || strcmpi(pOpt, 'all')
        pos = params.pos; % [Angstrom]
        spin0t = spin0(:,:,idx);
        spinTt = spinT(:,:,idx);
        spint = spin(:,:,idx);

        % move the orgin
        prompt = sprintf('Please select the origin (Dims: %d x %d x %d):\n', params.dims);
        orgn = input(prompt);
        if max(abs(orgn)) > abs(max(params.dims)); error('Origin out of range!'); end
        % convert to real unit (angstrom)
        pos = pos - orgn .* diag(ion.abc{ion.idx})';

        % define the spatial domain to be plotted
        prompt = sprintf('Please select the range (Dims: %d x %d x %d):\n', params.dims);
        rng = input(prompt);
        if rng < 0;  error('Range must be positive!'); end
        rng = rng .* diag(ion.abc{ion.idx})'; % convert to real unit (angstrom)
        inner = abs(pos(:,1)) <= rng(1) & abs(pos(:,2)) <= rng(2) & abs(pos(:,3)) <= rng(3);

        fig2 = figure;
        sf1 = subplot(1,3,1);
        quiver3(pos(inner,1), pos(inner,2), pos(inner,3), spin0t(inner,1), spin0t(inner,2), spin0t(inner,3), 'k');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Initial spin config at T = %.2f K, B = %.2f T', params.temp(tt), vecnorm(params.field(:,ff))))
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5 0.4 0.5 0.5]);

        sf2 = subplot(1,3,2);
        quiver3(pos(inner,1), pos(inner,2), pos(inner,3), spinTt(inner,1), spinTt(inner,2), spinTt(inner,3), 'b');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Thermalized spin config at T = %.2f K, B = %.2f T', params.temp(tt), vecnorm(params.field(:,ff))))

        sf3 = subplot(1,3,3);
        quiver3(pos(inner,1), pos(inner,2), pos(inner,3), spint(inner,1), spint(inner,2), spint(inner,3),'r');
        pbaspect([params.dims(1) params.dims(2) params.dims(3)])
        title(sprintf('Final spin config at T = %.2f K, B = %.2f T', params.temp(tt), vecnorm(params.field(:,ff))))
    end
    if strcmpi(pOpt, 'Mag_B') || strcmpi(pOpt, 'all')
        fg3 = figure;
        ax3 = fg3.CurrentAxes;
        hold on
        plot(vecnorm(params.field,1), squeeze(mean(meas.Mx(tt,:,:),3)), 'o', 'Color', colorOption{1});
        hold on
        plot(vecnorm(params.field,1), squeeze(mean(meas.My(tt,:,:),3)), 'o', 'Color', colorOption{2});
        plot(vecnorm(params.field,1), squeeze(mean(meas.Mz(tt,:,:),3)), 'o', 'Color', colorOption{3});
        legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
        xlabel('Magnetic Field (T)')
        ylabel('$\langle J \rangle$','Interpreter','latex')
        set(ax3, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'Mag_T') || strcmpi(pOpt, 'all')
        fg4 = figure;
        ax4 = fg4.CurrentAxes;
        hold on
        plot(params.temp, squeeze(mean(meas.Mx(:,ff,:),3)), 'o', 'Color', colorOption{1});
        hold on
        plot(params.temp, squeeze(mean(meas.My(:,ff,:),3)), 'o', 'Color', colorOption{2});
        plot(params.temp, squeeze(mean(meas.Mz(:,ff,:),3)), 'o', 'Color', colorOption{3});
        legend('$\langle J_x \rangle$','$\langle J_y \rangle$','$\langle J_z \rangle$','interpreter','latex');
        xlabel('Temperature (K)')
        ylabel('$\langle J \rangle$','Interpreter','latex')
        set(ax4, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'Cv_T') || strcmpi(pOpt, 'all')
        fg5 = figure;
        ax5 = fg5.CurrentAxes;
        hold on
        % E2m = squeeze(sum(meas.Em(:,ff).^2, 3)); % <E^2>
        % Em2 = squeeze(sum(meas.Em(:,ff), 3).^2); % <E>^2
        % plot(params.temp, (E2m - Em2)./ (kB * params.temp.^2), '-o', 'Color', colorOption{1});
        Em_var = squeeze(var(meas.Em,0,3))';
        Cv_T = Em_var./ (kB * params.temp.^2);
        plot(params.temp, Cv_T, '-o', 'Color', colorOption{1});
        ylabel('C$_v$(T)','Interpreter','latex')
        % plot(params.temp, log(Cv_T), '-o', 'Color', colorOption{1});
        % ylabel('log[C$_v$(T)]','Interpreter','latex')
        xlabel('Temperature (K)')
        set(ax5, 'FontSize', FntSz);
    end
    if strcmpi(pOpt, 'Cv_B') || strcmpi(pOpt, 'all')
        fg6 = figure;
        ax6 = fg6.CurrentAxes;
        hold on
        % E2m = squeeze(sum(meas.Em(tt,:).^2, 3)); % <E^2>
        % Em2 = squeeze(sum(meas.Em(tt,:), 3).^2); % <E>^2
        % plot(vecnorm(params.field,1), (E2m - Em2)./ (kB * vecnorm(params.field,1).^2), '-o', 'Color', colorOption{1});
        Em_var = squeeze(var(meas.Em,0,3));
        Cv_B = Em_var./ (kB * vecnorm(params.field,1).^2);
        plot(vecnorm(params.field,1), Cv_B, '-o', 'Color', colorOption{1});
        ylabel('C$_v$(B)','Interpreter','latex')
        % plot(vecnorm(params.field,1), log(Cv_B), '-o', 'Color', colorOption{1});
        % ylabel('log[C$_v$(B)]','Interpreter','latex')
        xlabel('Magnetic Field (T)')
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
                plt(tt, ff) = plot(1:cutoff, dEt(1:cutoff), 'LineWidth', 1.5);
                lgd{end+1} = ['Temp: ' num2str(params.temp(tt)) ', Field: ' num2str(ff)];
            end
        end
        xlabel('Iteration', 'Interpreter', 'latex');
        ylabel('$\Delta E$', 'Interpreter', 'latex');
        title('Energy Change Evolution over Iterations', 'Interpreter', 'latex');
        set(ax4, 'FontSize', 12);
        grid on;
        legend(lgd, 'Location', 'best');
    end
end
end