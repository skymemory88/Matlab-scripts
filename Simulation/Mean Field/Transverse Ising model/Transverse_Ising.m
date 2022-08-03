clear
format long g;

Options.plot = 1; % 1. 3D scatter plot; 2. 2D color map.

sigx = [0   1
        1   0];
sigy = [ 0   1i
        -1i  0];
sigz = [1   0
        0  -1];

Sx = 1/2.*sigx;
Sy = 1/2.*sigy;
Sz = 1/2.*sigz;

Ix = 1/2.*sigx;
Iy = 1/2.*sigy;
Iz = 1/2.*sigz;

% N = 100; % set interaction range
D = 1; % setting dimensions
k = 0.0862875; % Boltzmann coefficient
J = 0.6; % coupling strength

minField = 0;
maxField = 6;
hT = linspace(minField, maxField, 300);

minTemp = 0.01;
maxTemp = 3;
temp = linspace(minTemp, maxTemp, 150);

vector = zeros(length(hT)*length(temp),3);
map = zeros(length(hT), length(temp));

delta = 1e-5;
iterator = 1;
A = 0.0;
Zee = 0.08;
newAvSz = 0;

for ii = 1:length(hT)
    for jj = 1:length(temp)
        avSz = 0.1; %initial guess of average Sz
        beta = 1/(k*temp(jj));
        while true
            if A ~= 0 % Case: finite hyperfine interaction
%                 Hamlt = -J*avSz*kron(eye(2),Sz) + hT*kron(eye(2),Sx) + A.*(kron(Ix,Sx)+kron(Iy,Sy)+kron(Iz,Sz));
                Hamlt = -J*avSz*kron(eye(2),Sz) + Zee*hT(ii)*kron(eye(2),Sx) + A.*(kron(Ix,Sx)+kron(Iy,Sy)+kron(Iz,Sz));
                [wav, En] = eig(Hamlt);
                Z = trace(exp(-beta.*En)); % Calculate the partition function weight
                newAvSz = diag(exp(-beta*En))'*diag(wav'*kron(eye(2),Sz)*wav)/Z; % Direct product by two because of the two neighbours
            else % Case: no hyperfine interaction
                Hamlt = -J*avSz*Sz + Zee*hT(ii)*Sx;
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
        vector(iterator,:) = [temp(jj), hT(ii), avSz];
        map(ii,jj) = avSz;
        iterator = iterator + 1;
        %if mod(iterater, 1000) == 0;
        %    iterater, temp, hT
        %end
    end
end

Bc = zeros(length(temp),1); % critical fields (at fixed temperatre)
for ii = 1:length(temp)
    dif = diff(map(:,ii)); % crude differentiation
    [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', 2e-2);
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
    [~,idx] = findpeaks(-dif, 'Npeaks', 1, 'MinPeakProminence', 2e-2);
    if isempty(idx)
        Tc(ii) = 0;
    else
        Tc(ii) = temp(idx);
    end
end
BB = hT(Tc~=0)'; % remove zeros and nonsensical values
phB = [TT nonzeros(Bc); nonzeros(Tc) BB]; % Phase boundary

figure
box on
hold on
switch Options.plot
    case 1
        phase = scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
        view([0 90]);
    case 2
        temps = repmat(temp, length(hT), 1);
        fields = repmat(hT', 1, length(temp));
        phase = pcolor(temps, fields, map);
        phase.EdgeColor = 'none';
        plot(phB(:,1), phB(:,2), '.r');
end
set(gca, 'FontSize', 12)
xlabel('Temperature (K)')
ylabel('Magnetic Field (T)')