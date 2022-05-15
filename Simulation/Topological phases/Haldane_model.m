% Haldane model
t2min = 0;
t2max = 0.5;
mmin = 0;
mmax = 2*t2max;
lattc = 1.0; % lattice constant

fObj = figure;
ax = axes(fObj);
box on
hold(ax,'on');
[KX, KY, Er] = kBands(lattc, 0, 0);
uband = surf(ax, KX, KY, squeeze(Er(:,:,1)));
uband.EdgeColor = 'none';
lband = surf(ax,KX, KY, squeeze(Er(:,:,2)));
lband.EdgeColor = 'none';
xlabel(ax, '$\vec{k}_x$','Interpreter','latex')
ylabel(ax, '$\vec{k}_y$','Interpreter','latex')
zlabel(ax, 'Energy')
set(gca,'FontSize',14)
view(30,45)

slider1 = uicontrol('Parent', fObj, 'Style', 'slider',...
    'OuterPosition', [100 400 150 20], 'Value', 0, 'min', t2min, 'max', t2max);
sld1.label = uicontrol('Parent',fObj,'Style','text','Position', [250 400 40 20],...
                'String','Tnnn');
            
slider2 = uicontrol('Parent', fObj, 'Style', 'slider',...
    'OuterPosition', [100 380 150 20], 'Value', 0, 'min', mmin, 'max', mmax);
sld2.label = uicontrol('Parent',fObj,'Style','text', 'Position',[250 380 40 20],...
                'String','Mass');
            
slider1.Callback = @(es1,ed) band_update([uband lband], lattc, slider2.Value, es1.Value);
slider2.Callback = @(es2,ed) band_update([uband lband], lattc, es2.Value, slider1.Value);

function band_update(figs, lattc, mass, tnnn)
[~, ~, Er] = kBands(lattc, mass, tnnn);
figs(1).ZData = Er(:,:,1);
figs(2).ZData = Er(:,:,2);
drawnow
end

function [KX, KY, Er] = kBands(lattc, mass, tnnn)
tnn = 1.0; % nearest neighbour hopping amplitude

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
 
% k-vector
kx = linspace(-pi/lattc, pi/lattc, 200);
ky = linspace(-pi/lattc, pi/lattc, 200);
[KX, KY] = meshgrid(kx,ky);

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
        
        Hr = tnn*(sigx.*Hr1 - sigy.*Hr2) + sigz*(mass + Hr3); % total hamiltonian
        [~, ee] = eig(Hr);
        Er(ii,jj,:,1) = diag(ee);
    end
end
end