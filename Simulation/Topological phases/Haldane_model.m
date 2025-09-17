% Haldane model (Tnn defaults to 1)
mmin = 0.0;
mmax = 0.8;
t2min = -mmax;
t2max = +mmax;
phi = -pi/2; % hopping phase
lattc.Nk = 201; % number of points in k space

lattc.const = 1.0;
lattc.L = 25;
% unit vector
lattc.a = {[1 0].*lattc.const;
    [-1/2 -sqrt(3)/2].*lattc.const;
    [-1/2 sqrt(3)/2].*lattc.const};

lattc.b = {[0 -sqrt(3)].*lattc.const;
    [3/2 +sqrt(3)/2].*lattc.const;
    [-3/2 +sqrt(3)/2].*lattc.const};

lattc.a_zig = sqrt(3) * lattc.const; % zigzag period
lattc.a_arm = 3 * lattc.const; % armchair period

% k-vector
kx = linspace(-pi/lattc.const, pi/lattc.const, lattc.Nk);
ky = linspace(-pi/lattc.const, pi/lattc.const, lattc.Nk);
[KX, KY] = meshgrid(kx,ky);

% subfigure 1 (band diagram)
fobj1 = figure('Position', [50, 100, 1400, 440]);
ax1 = subplot('Position', [0.05 0.2 0.27 0.75]);
box on
hold(ax1,'on');
[Er, Hr] = kBands(kx, ky, lattc, 0, 0, phi);
uband = surf(ax1, KX, KY, squeeze(Er(:,:,1)).');
uband.EdgeColor = 'none';
lband = surf(ax1, KX, KY, squeeze(Er(:,:,2)).');
lband.EdgeColor = 'none';

xlabel(ax1, '\vec{k}_x','Interpreter','latex')
yticks([-pi -pi/2 0 pi/2 pi])
xticklabels(ax1, {'-\pi','-\pi/2','0','\pi/2','\pi'});

ylabel(ax1, '\vec{k}_y','Interpreter','latex')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels(ax1, {'-\pi','-\pi/2','0','\pi/2','\pi'});

zlabel(ax1, 'Energy')
set(ax1,'FontSize',15)
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
    [linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
    [ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows', 'stable');
cmap = [flip(cmap,1); cmap];
colormap(cmap)
view(30,35)

% subfigure 2 (ribbon plot-zigzag)
figure(fobj1);
ax2 = subplot('Position', [0.36 0.2 0.27 0.75]);
box on
hold(ax2,'on');

kx_zig = linspace(-pi/lattc.a_zig, pi/lattc.a_zig, lattc.Nk);   % BZ for zig-zag ribbon
En_zig = band_ribbon_zig(lattc, kx_zig, 0, 0, phi);

ribbon_zig = plot(ax2, kx_zig, En_zig, '-w');   % use black lines so they're visible

xline(ax2,  2*pi/(3*lattc.a_zig), ':');
xline(ax2, -2*pi/(3*lattc.a_zig), ':');
xlim(ax2, [-pi, pi]/lattc.a_zig);
xticks(ax2, (-pi:pi/2:pi)/lattc.a_zig);
xticklabels(ax2, {'-\pi/a_{\rm zig}','-\pi/(2a_{\rm zig})','0',...
    '\pi/(2a_{\rm zig})','\pi/a_{\rm zig}'});
xlabel(ax2, '\vec{k}_x','Interpreter','latex')
ylabel(ax2, 'E_n','Interpreter','latex')
set(ax2,'FontSize',15)

% subfigure 3 (ribbon plot-armchair)
figure(fobj1);
ax3 = subplot('Position', [0.67 0.2 0.27 0.75]);
box on;
hold(ax3, 'on');

kx_arm = linspace(-pi/lattc.a_arm, pi/lattc.a_arm, lattc.Nk);
En_arm = band_ribbon_armchair(lattc, kx_arm, 0, 0, phi);

ribbon_arm = plot(ax3, kx_arm, En_arm, '-w');                   % plot armchair ribbon bands

xlim(ax3, [-pi, pi]/lattc.a_arm);
xticks(ax3, (-pi:pi/2:pi)/lattc.a_arm);
xticklabels(ax3, {'-\pi/a_{\rm arm}','-\pi/(2a_{\rm arm})','0', ...
    '\pi/(2a_{\rm arm})','\pi/a_{\rm arm}'});
xlabel(ax3, '\vec{k}_x','Interpreter','latex');
ylabel(ax3, 'E_n','Interpreter','latex');
set(ax3,'FontSize',15);

% Slider setup (T2/T1 & Mass)
slider1 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'Position', [100 30 150 20], 'Value', 0.0, 'min', t2min, 'max', t2max);
sld1.label = uicontrol('Parent',fobj1,'Style','text','Position',  [310 30 40 20],...
    'String','T2/T1');
box1 = uicontrol('Style', 'edit','Position', [250 30 50 20], 'String', num2str(slider1.Value));

slider2 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'Position', [100 10 150 20], 'Value', 0.0, 'min', mmin, 'max', mmax);
sld2.label = uicontrol('Parent',fobj1,'Style','text', 'Position', [310 10 50 20],...
    'String','Mass/T1');
box2 = uicontrol('Style', 'edit','Position', [250 10 50 20], 'String', num2str(slider2.Value));

% slider update
slider1.Callback = @(es,ed) band_update({uband, lband, ribbon_zig, ribbon_arm}, kx, ky, lattc,...
    slider2.Value, es.Value, phi, box1, box2);
slider2.Callback = @(es,ed) band_update({uband, lband, ribbon_zig, ribbon_arm}, kx, ky, lattc,...
    es.Value, slider1.Value, phi, box1, box2);

% helper function for plot live update
function band_update(figs, kx, ky, lattc, mass, tnnn, phi, box1, box2)
set(box1, 'String', num2str(tnnn));
set(box2, 'String', num2str(mass));

kx_zig = linspace(-pi/lattc.a_zig, pi/lattc.a_zig, lattc.Nk); % BZ for zig-zag ribbon
kx_arm = linspace(-pi/lattc.a_arm, pi/lattc.a_arm, lattc.Nk); % BZ for armchair ribbon

[Er, ~] = kBands(kx, ky, lattc, mass, tnnn, phi);
En_zig = band_ribbon_zig(lattc, kx_zig, mass, tnnn, phi);
En_arm = band_ribbon_armchair(lattc, kx_arm, mass, tnnn, phi);

% update the band diagram
figs{1}.ZData = squeeze(Er(:,:,1)).';
figs{2}.ZData = squeeze(Er(:,:,2)).';

% update the ribbon plots
for ii = 1:numel(figs{3})
    figs{3}(ii).YData = real(En_zig(:,ii));
end
for ii = 1:numel(figs{4})
    figs{4}(ii).YData = real(En_arm(:,ii));
end
drawnow
end

function [Er, Hr] = kBands(kx, ky, lattc, mass, tnnn, phi)
% Pauli matrices
sigx = [0   1;  1   0];
sigy = [0 -1i;  1i  0];
sigz = [1   0;  0  -1];

% unit vector
a = lattc.a;
b = lattc.b;

Hr = double.empty(length(kx), length(ky), 2, 2, 0);
Er = double.empty(length(kx), length(ky), 2, 0);
for ii = 1:length(kx)
    for jj = 1:length(ky)
        kv = [kx(ii) ky(jj)]; % k-vector

        Hr1 = 0; Hr2 = 0;
        for kk = 1:length(a) % N.N. contribution
            Hr1 = Hr1 + cos(kv*a{kk}');
            Hr2 = Hr2 + sin(kv*a{kk}');
        end
        Hr3 = 0; Hr0 = 0;
        for kk = 1:length(b) % N.N.N. contribution
            Hr0 = Hr0 + 2*tnnn*cos(phi)*cos(kv*b{kk}');
            Hr3 = Hr3 - 2*tnnn*sin(phi)*sin(kv*b{kk}');
        end

        Hr(ii,jj,:,:,1) = (Hr1.*sigx - Hr2.*sigy) + (mass + Hr3).*sigz...
            + Hr0.*eye(2); % total hamiltonian
        [~, ee] = eig(squeeze(Hr(ii,jj,:,:,1)));
        Er(ii,jj,:,1) = diag(ee);
    end
end
end

function En = band_ribbon_zig(lattc, kx, mass, tnnn, phi)
% Zigzag ribbon (periodic along x with a_zig = sqrt(3)*a, open along y).
% Basis: [A1,B1,A2,B2,...,AN,BN], N = lattc.L.
% Couplings:
%   A_j <-> B_j          :  2 t1 cos(k a_zig / 2)
%   B_j <-> A_{j+1}      :  t1
%   A_j <-> A_{j+1}      :  2 t2 cos(k a_zig / 2 + phi)
%   B_j <-> B_{j+1}      :  2 t2 cos(k a_zig / 2 - phi)
%   Onsite: A = -m + 2 t2 cos(k a_zig - phi),  B = +m + 2 t2 cos(k a_zig + phi)

tnn   = 1.0;
N    = lattc.L;
D    = 2*N;

En = zeros(numel(kx), D);
for p = 1:numel(kx)
    k = kx(p);
    H = zeros(D, D);

    % On-site (A/B)
    for j = 1:N
        iA = 2*j - 1; iB = 2*j;
        H(iA,iA) = -mass + 2*tnnn*cos(k*lattc.a_zig - phi);
        H(iB,iB) =  mass + 2*tnnn*cos(k*lattc.a_zig + phi);
    end

    % Nearest neighbors
    c_k = 2*tnn*cos(k*lattc.a_zig/2);
    for j = 1:N
        iA = 2*j - 1; iB = 2*j;
        % A_j <-> B_j  (k-dependent)
        H(iA,iB) = H(iA,iB) + c_k; H(iB,iA) = H(iA,iB);

        % B_j <-> A_{j+1}
        if j < N
            iA2 = 2*(j+1) - 1;
            H(iB,iA2) = H(iB,iA2) + tnn; H(iA2,iB) = H(iB,iA2);
        end
    end

    % Next-nearest neighbors (same sublattice)
    cA = 2*tnnn*cos(k*lattc.a_zig/2 + phi);
    cB = 2*tnnn*cos(k*lattc.a_zig/2 - phi);
    for j = 1:N-1
        iA = 2*j - 1; iB = 2*j;
        iA2 = 2*(j+1) - 1; iB2 = 2*(j+1);

        % A_j <-> A_{j+1}
        H(iA,iA2) = H(iA,iA2) + cA; H(iA2,iA) = H(iA,iA2);
        % B_j <-> B_{j+1}
        H(iB,iB2) = H(iB,iB2) + cB; H(iB2,iB) = H(iB,iB2);
    end

    H = (H + H')/2;              % numerical hermiticity
    ev = eig(H);
    En(p,:) = sort(real(ev)).';
end
end

function En_arm = band_ribbon_armchair(lattc, kx_arm, mass, tnnn, phi)
% Armchair ribbon (periodic along x with a_arm = 3*a, open along y).
% Basis: [A1, B1, A2, B2, ..., A_N, B_N], N = lattc.L (number of dimer rows).
% Couplings:
%   A_j <-> B_j (same row, armchair direction):  2 t1 cos(k_x * a_arm / 2)
%   A_j <-> B_{j+1} (downward neighbor):  t1    (direct vertical/slanted bond)
%   B_j <-> A_{j+1} (downward neighbor):  t1    (direct vertical/slanted bond)
%   On-site: A_j = -m + 2 t2 cos(k_x * a_arm - phi),
%            B_j = +m + 2 t2 cos(k_x * a_arm + phi)
%   A_j <-> A_{j+1} (next-nearest between rows):  2 t2 cos(k_x * a_arm/2 + phi)
%   B_j <-> B_{j+1} (next-nearest between rows):  2 t2 cos(k_x * a_arm/2 - phi)

tnn = 1.0;
N  = lattc.L;
D  = 2 * N;
En_arm = zeros(numel(kx_arm), D);
a_arm = 3 * lattc.const;
for p = 1:numel(kx_arm)
    k = kx_arm(p);
    H = zeros(D, D);
    % On-site energies (including mass and t2 flux terms)
    for j = 1:N
        iA = 2*j - 1; iB = 2*j;
        H(iA,iA) = -mass + 2*tnnn*cos(k * a_arm - phi);
        H(iB,iB) =  mass + 2*tnnn*cos(k * a_arm + phi);
    end
    % Nearest-neighbor hoppings (t1)
    c_arm = 2 * tnn * cos(k * a_arm/2);   % armchair horizontal A–B coupling
    for j = 1:N
        iA = 2*j - 1; iB = 2*j;
        % A_j <-> B_j (same-row armchair bond with k-phase)
        H(iA,iB) = H(iA,iB) + c_arm;
        H(iB,iA) = H(iA,iB);  % symmetric
        % Inter-row neighbors:
        if j < N
            % B_j <-> A_{j+1}
            iA2 = 2*(j+1) - 1;
            H(iB,iA2) = H(iB,iA2) + tnn;
            H(iA2,iB) = H(iB,iA2);
            % A_j <-> B_{j+1}
            iB2 = 2*(j+1);
            H(iA,iB2) = H(iA,iB2) + tnn;
            H(iB2,iA) = H(iA,iB2);
        end
    end
    % Next-nearest neighbors (t2, same sublattice couplings)
    cA = 2 * tnnn * cos(k * a_arm/2 + phi);
    cB = 2 * tnnn * cos(k * a_arm/2 - phi);
    for j = 1:N-1
        iA = 2*j - 1;   iB = 2*j;
        iA2 = 2*(j+1) - 1; iB2 = 2*(j+1);
        % A_j <-> A_{j+1}
        H(iA, iA2) = H(iA, iA2) + cA;
        H(iA2, iA) = H(iA, iA2);
        % B_j <-> B_{j+1}
        H(iB, iB2) = H(iB, iB2) + cB;
        H(iB2, iB) = H(iB, iB2);
    end
    H = (H + H')/2;            % ensure Hermiticity
    ev = eig(H);
    En_arm(p, :) = sort(real(ev)).';
end
end