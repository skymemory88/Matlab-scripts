% Haldane model (Tnn defaults to 1)
mmin = 0.0;
mmax = 0.5;
t2min = -mmax;
t2max = +mmax;
phi = -pi/2; % hopping phase

lattc.const = 1.0;
lattc.L = 25;
% unit vector
lattc.a = {[1 0].*lattc.const;
    [-1/2 -sqrt(3)/2].*lattc.const;
    [-1/2 sqrt(3)/2].*lattc.const};

lattc.b = {[0 -sqrt(3)].*lattc.const;
     [3/2 +sqrt(3)/2].*lattc.const;
     [-3/2 +sqrt(3)/2].*lattc.const};

% k-vector
kx = linspace(-pi/lattc.const, pi/lattc.const, 200);
ky = linspace(-pi/lattc.const, pi/lattc.const, 200);
[KX, KY] = meshgrid(kx,ky);

fobj1 = figure('Position', [50, 100, 1200, 440]);
ax1 = subplot('Position', [0.1 0.2 0.35 0.75]);

box on
hold(ax1,'on');
[Er, Hr] = kBands(kx, ky, lattc, 0, 0, phi);
uband = surf(ax1, KX, KY, squeeze(Er(:,:,1)));
uband.EdgeColor = 'none';
lband = surf(ax1,KX, KY, squeeze(Er(:,:,2)));
lband.EdgeColor = 'none';
xlabel(ax1, '$\vec{k}_x$','Interpreter','latex')
xticks([-pi -pi/2 0 pi/2 pi])
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
ylabel(ax1, '$\vec{k}_y$','Interpreter','latex')
yticks([-pi -pi/2 0 pi/2 pi])
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
zlabel(ax1, 'Energy')
set(ax1,'FontSize',15)
cmap = unique([[0 0 0];[zeros(20,1),linspace(0,0.5,20)',linspace(0.4,1,20)'];...
[linspace(0,1,36)',0.5*ones(36,1),linspace(1,0,36)'];...
[ones(30,1),linspace(0.5,1,30)',linspace(0,1,30)']],'rows');
cmap = [flip(cmap,1); cmap];
colormap(cmap)
view(30,35)

figure(fobj1);
ax2 = subplot('Position', [0.6 0.2 0.35 0.75]);
box on
hold(ax2,'on');
En = band_ribbon(lattc, kx, squeeze(Hr(100,:,:,:)), 0, phi);
cband = plot(KX(1,:), squeeze(En(:,:)), '-k');
xlabel(ax2, '$\vec{k}_x$','Interpreter','latex')
ylabel(ax2, '$E_n$','Interpreter','latex')
set(ax2,'FontSize',15)

slider1 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'OuterPosition', [100 30 150 20], 'Value', 0.0, 'min', t2min, 'max', t2max);
sld1.label = uicontrol('Parent',fobj1,'Style','text','Position',  [310 30 40 20],...
                'String','T2/T1');
box1 = uicontrol('Style', 'edit','Position', [250 30 50 20], 'String', num2str(slider1.Value));
            
slider2 = uicontrol('Parent', fobj1, 'Style', 'slider',...
    'OuterPosition', [100 10 150 20], 'Value', 0.0, 'min', mmin, 'max', mmax);
sld2.label = uicontrol('Parent',fobj1,'Style','text', 'Position', [310 10 50 20],...
                'String','Mass/T1');
box2 = uicontrol('Style', 'edit','Position', [250 10 50 20], 'String', num2str(slider2.Value));
            
slider1.Callback = @(es,ed) band_update({uband lband cband}, kx, ky, lattc,...
    slider2.Value, es.Value, phi, box1, box2);
% box1.Callback = @(es,ed) band_update({uband lband cband}, kx, ky, lattc,...
%     slider2.Value, str2double(es.Value), phi, chainL, box1, box2);

slider2.Callback = @(es,ed) band_update({uband lband cband}, kx, ky, lattc,...
    es.Value, slider1.Value, phi, box1, box2);
% box2.Callback = @(es,ed) band_update({uband lband cband}, kx, ky, lattc,...
%     str2double(es.Value), slider1.Value, phi, chainL, box1, box2);

function band_update(figs, kx, ky, lattc, mass, tnnn, phi, box1, box2)
set(box1, 'String', num2str(tnnn));
set(box2, 'String', num2str(mass));
[Er, Hr] = kBands(kx, ky, lattc, mass, tnnn, phi);
En = band_ribbon(lattc, kx, squeeze(Hr(100,:,:,:)), tnnn, phi);
figs{1}.ZData = Er(:,:,1);
figs{2}.ZData = Er(:,:,2);
for ii = 1:size(figs{3},1)
    figs{3}(ii).YData = real(squeeze(En(:,ii)));
end
drawnow
end

function En = band_ribbon(lattc, kv, Hr, tnnn, phi)
% unit vector
a = lattc.a;
b = lattc.b;
NN = lattc.L;

En = double.empty(size(Hr,1), NN*2, 0);
for ii = 1:size(Hr,1)
    hr = squeeze(Hr(ii,:,:));
    ham = kron(eye(NN), hr/2);
    
% NN hopping betwen cells
    ham((NN+1)*2+1 : (2*NN+1)*4 : end) = cos([kv(ii) 0]*a{2}') + 1i*sin([kv(ii) 0]*a{2}'); % NN hopping betwen cells
    ham((NN+1)*2+(2*NN+1)*2+1 : (2*NN+1)*4 : end) = cos([kv(ii) 0]*a{3}') + 1i*sin([kv(ii) 0]*a{3}'); % NN hopping betwen cells
    
    ham(3 : 4*(2*NN+1) : end-NN*2) = tnnn*exp(-1i*phi)*(cos([kv(ii) 0]*b{1}') + 1i*sin([kv(ii) 0]*b{1}')); % NNN hopping betwen cells
    ham(3+2*NN+1 : 4*(2*NN+1) : end) = tnnn*exp(-1i*phi)*(cos([kv(ii) 0]*b{1}') + 1i*sin([kv(ii) 0]*b{1}')); % NNN hopping betwen cells
    ham(3+2*(2*NN+1) : 4*(2*NN+1) : end) = tnnn*exp(-1i*phi)*(cos([kv(ii) 0]*b{3}') + 1i*sin([kv(ii) 0]*b{3}')); % NNN hopping betwen cells
    ham(3+3*(2*NN+1) : 4*(2*NN+1) : end) = tnnn*exp(-1i*phi)*(cos([kv(ii) 0]*b{3}') + 1i*sin([kv(ii) 0]*b{3}')); % NNN hopping betwen cells

    ham = ham' + ham; % alternative way of constructing the hamiltonian
    ham(1) = ham(1)-tnnn*sin(phi)*sin(phi)*sin([kv(ii) 0]*a{2}')...
        - tnnn*sin(phi)*sin([kv(ii) 0]*a{3}'); % account for missing bonds at the edges
    ham(end) = ham(end)-tnnn*sin(phi)*sin(phi)*sin([kv(ii) 0]*a{2}')...
        - tnnn*sin(phi)*sin([kv(ii) 0]*a{3}'); % account for missing bonds at the edges
    [~, ee] = eig(ham);
    En(ii,:,1) = real(diag(ee));
end
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