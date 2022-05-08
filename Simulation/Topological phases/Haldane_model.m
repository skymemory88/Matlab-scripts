% Haldane model
lattc = 1; % lattice constant
tnn = 1.0; % nearest neighbour hopping amplitude
tnnn = 0.0; % next nearest neighbour hopping amplitude
mass = 0.2; % mass term

% unit vector
a = {[1 0].*lattc;
    [-1/2 -sqrt(3)/2].*lattc;
    [-1/2 sqrt(3)/2].*lattc};

b = {[3/2 -sqrt(3)/2].*lattc;
     [3/2 +sqrt(3)/2].*lattc;
     [-3/2 +sqrt(3)/2].*lattc;
     [-3/2 -sqrt(3)/2].*lattc
     [0 -sqrt(3)].*lattc
     [0 +sqrt(3)].*lattc};

% k-vector
kx = linspace(-pi/lattc, pi/lattc, 200);
ky = linspace(-pi/lattc, pi/lattc, 200);
[KX, KY] = meshgrid(kx,ky);

% Pauli matrices
sigx = [0 1; 1 0];
sigy = [0 -1i; 1i 0];
sigz = [1 0; 0 -1];

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
            Hr3 = Hr3 + sin(kv*b{kk}');
        end
        
        Hr = tnn*(sigx.*Hr1 - sigy.*Hr2) + sigz*(mass + tnnn*Hr3); % total hamiltonian
        [~, ee] = eig(Hr);
        Er(ii,jj,:,1) = diag(ee);
    end
end

figure
box on
hold on
uband = surf(KX, KY, squeeze(Er(:,:,1)));
uband.EdgeColor = 'none';
lband = surf(KX, KY, squeeze(Er(:,:,2)));
lband.EdgeColor = 'none';
xlabel('$\vec{k}_x$','Interpreter','latex')
ylabel('$\vec{k}_y$','Interpreter','latex')
zlabel('Energy')
set(gca,'FontSize',14)
view(30,45)