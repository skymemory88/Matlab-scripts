% Haldane model (Tnn defaults to 1)
mmin = 0;
mmax = 0.5;
t2min = -0.5*mmax;
t2max = +0.5*mmax;
lattc = 1.0; % lattice constant
chainL = 20;
 
% k-vector
kx = linspace(-pi/lattc, pi/lattc, 200);
ky = linspace(-pi/lattc, pi/lattc, 200);
[KX, KY] = meshgrid(kx,ky);

fobj1 = figure('Position', [50, 100, 1200, 440]);
ax1 = subplot('Position', [0.1 0.2 0.35 0.75]);

box on
hold(ax1,'on');
[Er, Hr] = kBands(kx, ky, lattc, 0, 0);
uband = surf(ax1, KX, KY, squeeze(Er(:,:,1)));
uband.EdgeColor = 'none';
lband = surf(ax1,KX, KY, squeeze(Er(:,:,2)));
lband.EdgeColor = 'none';
xlabel(ax1, '$\vec{k}_x$','Interpreter','latex')
ylabel(ax1, '$\vec{k}_y$','Interpreter','latex')
zlabel(ax1, 'Energy')
set(ax1,'FontSize',12)
view(30,35)

figure(fobj1);
ax2 = subplot('Position', [0.6 0.2 0.35 0.75]);
box on
hold(ax2,'on');
En = ribbon(squeeze(Hr(100,:,:,:)), chainL);
cband = plot(KY(:,1), squeeze(En(:,:)), '-k');
xlabel(ax1, '$\vec{k}_y$','Interpreter','latex')
ylabel(ax1, '$E_n$','Interpreter','latex')
set(ax2,'FontSize',12)

slider1 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'OuterPosition', [100 400 150 20], 'Value', 0, 'min', t2min, 'max', t2max);
sld1.label = uicontrol('Parent',fobj1,'Style','text','Position', [250 400 40 20],...
                'String','T2/T1');
            
slider2 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'OuterPosition', [100 380 150 20], 'Value', 0.25, 'min', mmin, 'max', mmax);
sld2.label = uicontrol('Parent',fobj1,'Style','text', 'Position',[250 380 40 20],...
                'String','Mass/T1');
            
slider1.Callback = @(es1,ed) band_update({uband lband cband}, kx, ky, lattc,...
    slider2.Value, es1.Value, chainL);
slider2.Callback = @(es2,ed) band_update({uband lband cband}, kx, ky, lattc,...
    es2.Value, slider1.Value, chainL);

function band_update(figs, kx, ky, lattc, mass, tnnn, NN)
[Er, Hr] = kBands(kx, ky, lattc, mass, tnnn);
En = ribbon(squeeze(Hr(100,:,:,:)), NN);
figs{1}.ZData = Er(:,:,1);
figs{2}.ZData = Er(:,:,2);
for ii = 1:size(figs{3},1)
    figs{3}(ii).YData = squeeze(En(:,ii));
end
drawnow
end

function En = ribbon(Hr, NN)
En = double.empty(size(Hr,1), NN*2, 0);
for ii = 1:size(Hr,1)
    hr = squeeze(Hr(ii,:,:));
    ham = kron(eye(NN),hr);
    ham(NN*2+3:NN*4+2:end) = 2;
    ham = (ham'+ham)/2;
    [~, ee] = eig(ham);
    En(ii,:,1) = diag(ee);
end
end

function [Er, Hr] = kBands(kx, ky, lattc, mass, tnnn)
% Pauli matrices
sigx = [0   1;  1   0];
sigy = [0 -1i;  1i  0];
sigz = [1   0;  0  -1];

% unit vector
a = {[1 0].*lattc;
    [-1/2 -sqrt(3)/2].*lattc;
    [-1/2 sqrt(3)/2].*lattc};

b = {[0 -sqrt(3)].*lattc;
     [3/2 +sqrt(3)/2].*lattc;
     [-3/2 +sqrt(3)/2].*lattc};

Hr = double.empty(length(kx), length(ky), 2, 2, 0);
Er = double.empty(length(kx), length(ky), 2, 0);
for ii = 1:length(kx)
    for jj = 1:length(ky)
        kv = [kx(ii) ky(jj)]; % k-vector
        
        Hr1 = 0; Hr2 = 0; Hr3 = 0;
        for kk = 1:length(a) % N.N. contribution
            Hr1 = Hr1 + cos(kv*a{kk}');
            Hr2 = Hr2 + sin(kv*a{kk}');
        end
        
        for kk = 1:length(b) % N.N.N. contribution
            Hr3 = Hr3 + 2*tnnn*sin(kv*b{kk}');
        end
        
        Hr(ii,jj,:,:,1) = (sigx.*Hr1 - sigy.*Hr2) + sigz*(mass + Hr3); % total hamiltonian
        [~, ee] = eig(squeeze(Hr(ii,jj,:,:,1)));
        Er(ii,jj,:,1) = diag(ee);
    end
end
end