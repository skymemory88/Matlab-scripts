clearvars
format long g;

Options.plot = 2; % 1. 3D scatter plot; 2. 2D color map.
Option.pProm = 1e-2; % minimum peak prominance for phase boundary search

sigx = [0   1
        1   0];
sigy = [ 0   1i
        -1i  0];
sigz = [1   0
        0  -1];

J = 1/2; % electronic spin length
Sx = J * sigx;
Sy = J * sigy;
Sz = J * sigz;

I = 1/2; % nuclear spin length
Ix = I * sigx;
Iy = I * sigy;
Iz = I * sigz;

% N = 100; % set interaction range
D = 3; % setting dimensions
kB = 0.0862875; % Boltzmann coefficient
% Jij = -0.585; % coupling strength
Jij = -1; % coupling strength
eMeas = abs(Jij); % energy normalize agains coupling strength

minField = 0.00;
maxField = 0.60;
hT = linspace(minField, maxField, 600);

minTemp = 0.01;
maxTemp = 3.0;
temp = linspace(minTemp, maxTemp, 800);

vector = zeros(length(hT)*length(temp),3);
map = zeros(length(hT), length(temp));

delta = 1e-5;
A = 0.0;
% Zee = 0.08; % scale the Zeeman interaction strength
Zee = 1;
newAvSz = 0;

for ii = 1:length(hT)
    parfor jj = 1:length(temp)
        avSz = 0.1; %initial guess of average Sz
        beta = 1/(kB*temp(jj));
        while true
            if A ~= 0 % Case: finite hyperfine interaction
%                 Hamlt = Jij*avSz*kron(eye(2),Sz) - hT*kron(eye(2),Sx) + A.*(kron(Ix,Sx)+kron(Iy,Sy)+kron(Iz,Sz));
                Hamlt = Jij*avSz*kron(eye(2),Sz) - Zee*hT(ii)*kron(eye(2),Sx) + A.*(kron(Ix,Sx)+kron(Iy,Sy)+kron(Iz,Sz));
                [wav, En] = eig(Hamlt);
                Z = trace(exp(-beta.*En)); % Calculate the partition function weight
                newAvSz = diag(exp(-beta*En))'*diag(wav'*kron(eye(2),Sz)*wav)/Z; % Direct product by two because of the two neighbours
            else % Case: no hyperfine interaction
                Hamlt = Jij*avSz*Sz - Zee*hT(ii)*Sx;
                [wav, En] = eig(Hamlt);
                Z = trace(exp(-beta.*En)); % Calculate the partition function weight
                newAvSz = diag(exp(-beta*En))'*diag(wav'*Sz*wav)/Z; % Calculate the expectation value of Jz
            end
            if abs(newAvSz - avSz) < delta
                break
            else
                avSz = newAvSz;
            end
        end
        vector(jj,:) = [temp(jj), hT(ii), avSz];
        map(ii,jj) = avSz;
    end
end

Bc = zeros(length(temp),1); % critical fields (at fixed temperatre)
for ii = 1:length(temp)
    dif = diff(map(:,ii)); % crude differentiation
    [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Option.pProm);
    if isempty(idx)
        Bc(ii) = 0;
    else
        Bc(ii) = hT(idx);
    end
end
TT = temp(Bc~=0)'; % remove zeros and nonsensical values

Tc = zeros(length(hT),1); % critical temperatures (at fixed field)
for ii = 1:length(hT)
    dif = diff(map(ii,:)); % crude differentiation
    [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', Option.pProm);
    if isempty(idx)
        Tc(ii) = 0;
    else
        Tc(ii) = temp(idx);
    end
end
BB = hT(Tc~=0)'; % remove zeros and nonsensical values
phB = [kB*TT nonzeros(Bc); kB*nonzeros(Tc) BB]./eMeas; % Phase boundary

figure
box on
hold on
switch Options.plot
    case 1
        phase = scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
        view([0 90]);
        plot(phB(:,1), phB(:,2), zeros(size(phB,1),1), '.r');
    case 2
        temps = repmat(kB*temp./eMeas, length(hT), 1);
        fields = repmat(hT'./eMeas, 1, length(temp));
        phase = pcolor(temps, fields, map);
        phase.EdgeColor = 'none';
        plot(phB(:,1), phB(:,2), '.r');
end
set(gca, 'FontSize', 12)
xlim([0 maxTemp/eMeas])
ylim([0 maxField/eMeas])
xlabel('|T/J_{ij}|')
ylabel('|B/J_{ij}|')