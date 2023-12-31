Temperatures = 0.18;
theta = 0.0;
phi = -13.5;
nZee = true;

N_level = 8; % upper energy level to include in the calculation
gN = 1.192; % gN = mu/mu_N/<I> = 4.173/(7/2)
gLande = 1.25;
ionJ = 8; % Ho
ionI = 3.5; % Ho
lattice = [5.175 0 0;
           0 5.175 0;
           0 0 10.75]; % Lattice constant for LiHoF4
freq_total = linspace(3.54,3.76,200);
hbar = 1.055E-34; % Reduced Planck constant [J.s]
muN = 5.05078e-27; % Nuclear magneton [J/T]
J2meV = 6.24151e+21; % Joule to meV
f2E = hbar*2*pi*10^9*J2meV;% GHz to meV
muB = 9.274e-24; %[J/T]
mu0 = 4e-7*pi; % [H/m]
ELEf = gLande * muB; % Lande factor * Bohr magneton (J/T)
NUCf = gN * muN; % (J/T)

f_cav = 3.645; % cavity frequency (GHz)
filFctr = 0.0127*2.1;
rho = 4e30/det(lattice); % Holmium (magnetic moment) number density [m^-3]
gw0 = sqrt(mu0 * 2*pi * f_cav*1e9 * rho/2) * filFctr; % susceptibility prefactor [T^2/J. rad/s]^1/2
gw2 = gw0^2 * 2*pi * 1e-9; % [T^2/J. GHz]

if nZee == true
    nZee_path = 'Hz_I=1';
else
    nZee_path = 'Hz_I=0';
end

if strcmp(pathsep, ':')
    platform = 'Unix';
else
    platform = 'Win';
end

switch platform
    case 'Win'
        location = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\Matlab',...
            '\Susceptibilities\',nZee_path];
    case 'Unix'
        Options.location = ['/Users/yikaiyang/Library/CloudStorage/GoogleDrive-yikai.yang@epfl.ch/My Drive/'...
            'File sharing/PhD program/Research projects/LiHoF4 project/Data/',...
            'Simulations/MATLAB/Susceptibilities/', nZee_path];
end

filename = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=1.00', Temperatures, theta, phi),'.mat');
file = fullfile(location, filename);
load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
fields = vecnorm(fff,2,1);
eigenE = eee - min(eee,[],2); % normalize the eigen-energies
eigenW = vvv; % eigen-functions
temperature = repmat(ttt,length(fields),1);

%Initiate ionJ operators
Jz=diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
JhT.z=kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
Jmh=kron(Jm,eye(2*ionI+1));
JhT.x=(Jph+Jmh)/2;
JhT.y=(Jph-Jmh)/2i;

%Initiate I operators
Iz=diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
IhT.z=kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
Imh=kron(eye(2*ionJ+1),Im);
IhT.x=(Iph+Imh)/2;
IhT.y=(Iph-Imh)/2i;

% Single out <Jz+Iz> calculations
tz = double.empty(size(eigenW,2),size(eigenW,2),size(eigenW,1),0); % Expectation value of J pseudo-spin
ttz = double.empty(size(eigenW,2),size(eigenW,2),size(eigenW,1),0); % Expectation value of J-I pseudo-spin
Gc2 = double.empty(size(eigenW,2),size(eigenW,2),size(eigenW,1),0); % Coupling strength
for kk = 1:size(fields,2) % calculate susceptibility for all fields
    en = squeeze(eigenE(kk,:)); % Obtain the corresponding eigen energies [meV]
    if temperature(kk) ~= 0
        beta = 11.6/temperature(kk); %[meV^-1]
        zn = sum(exp(-beta*en));
        Z = exp(-beta*en)/zn;
        [n,np] = meshgrid(Z,Z);
        NN = n-np;
    else
        zn = zeros(size(en));
        zn(1) = 1;
        [n,np] = meshgrid(zn,zn);
        NN = n-np;
    end

    v = squeeze(eigenW(kk,:,:)); % Obtain the corresponding eigen vectors
    JzhT = JhT.z * ELEf;
    IzhT = IhT.z * NUCf;
    tz(:,:,kk,1) = v'  * JzhT * v;
    ttz(:,:,kk,1)  = v'  * (JzhT+IzhT) * v;
    Gc2(:,:,kk,1) = gw2 * ttz(:,:,kk,1) .* (ttz(:,:,kk,1).') .* NN;
end

figure;
hold on
for ii = 1:N_level-1
    plot(fields, sqrt(squeeze(abs(Gc2(ii,ii+1,:)))),'.');
end
xlabel('Magnetic Field (T)')
ylabel('Coupling strength (GHz)')