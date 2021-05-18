% function [temp,intc,Delt,eni,enj,J0i,J0j,Jsj,Jsi] = magn_magn_cpl(tempf,intf,scale)
% Toy model: Two coupled/uncoupled two-level-system (TLS)/multi-level-system (MLS) from spin-spin interaction
if ~exist('tempf','var')
    tempf = 0.0;
end
if ~exist('intf','var')
    intf = 0.5;
end
if ~exist('scale','var')
    scale = 2.0;
end
clearvars -except tempf intf scale
Delt = -1; % Ferromagnetic gap
intc = intf*Delt; % Spin-spin interaction strength
temp = -tempf*Delt; % Temperature
beta = 1/temp; % 1/kBT
Hx = linspace(0,-5*Delt,80);

% J = 8; % Ho
% J = 1; % Ho
J = 1/2;
Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jx=(Jp+Jm)/2;
% Jy=(Jp-Jm)/2i;

Jzhi = kron(eye(2*J+1),Jz); % Expand the Hilbert space (TLS-1)
Jphi = kron(eye(2*J+1),Jp);
Jmhi = kron(eye(2*J+1),Jm);
Jxhi = (Jphi+Jmhi)/2;
Jyhi = (Jphi-Jmhi)/2i;

Jzhj = kron(Jz,eye(2*J+1)); % Expand the Hilbert space (TLS-2)
Jphj = kron(Jp,eye(2*J+1));
Jmhj = kron(Jm,eye(2*J+1));
Jxhj = (Jphj+Jmhj)/2;
Jyhj = (Jphj-Jmhj)/2i;

% containers for two-level-system 1
eni = zeros(length(Jz),length(Hx)); % Eigen-energy
wavi = zeros(length(Jz),length(Jz),length(Hx)); % Eigen-state
J0i = zeros(length(Jz),length(Hx)); % Spin(z) expectation
% containers for two-level-system 2
enj = zeros(length(Jz),length(Hx));
wavj = zeros(length(Jz),length(Jz),length(Hx));
J0j = zeros(length(Jz),length(Hx));
% containers for the coupled two-level-systems
ens = zeros(length(Jzhi),length(Hx));
wav = zeros(length(Jzhi),length(Jzhi),length(Hx));
Jsi = zeros(length(Jzhi),length(Hx));
Jsj = zeros(length(Jzhj),length(Hx));
for ii = 1:length(Hx)
    Ham0i = scale*Delt*Jz + Hx(ii)*Jx;
    [v,e] = eig(Ham0i); % Diagonalize the hamiltonian
    e = real(diag(e)); % Take only the real part of the eigen-energy to form a diaganol matrix
    e = e-min(e); % Normalize the energy amplitude to the lowest eigen-energy
    [eni(:,ii),n] = sort(e); % sort the energy from lowest to the highest
    wavi(:,:,ii) = v(:,n); % sort the eigen-vectors in its basis accordingly
    zi = exp(-e*beta)/sum(exp(-e*beta)); % Partition function
    

    Ham0j = Delt/scale*Jz + Hx(ii)*Jx;
    [v,e] = eig(Ham0j); % Diagonalize the hamiltonian
    e = real(diag(e)); % Take only the real part of the eigen-energy to form a diaganol matrix
    e = e-min(e); % Normalize the energy amplitude to the lowest eigen-energy
    [enj(:,ii),n] = sort(e); % sort the energy from lowest to the highest
    wavj(:,:,ii) = v(:,n); % sort the eigen-vectors in its basis accordingly
    zj = exp(-e*beta)/sum(exp(-e*beta)); % Partition function
    
    
    Hami = scale*Delt*Jzhi + Hx(ii)*Jxhi;
    Hamj = Delt/scale*Jzhj + Hx(ii)*Jxhj;
    Hint = intc*Jzhi*Jzhj; % Ising interaction
%     Hint = intc*(Jxhi*Jxhj+Jyhi*Jyhj+Jzhi*Jzhj); % Heisenberg interaction
    Ham = Hami + Hamj + Hint;
    [v,e] = eig(Ham); % Diagonalize the hamiltonian
    e = real(diag(e)); % Take only the real part of the eigen-energy to form a diaganol matrix
    e = e-min(e); % Normalize the energy amplitude to the lowest eigen-energy
    [ens(:,ii),n] = sort(e); % sort the energy from lowest to the highest
    wav(:,:,ii) = v(:,n); % sort the eigen-vectors in its basis accordingly
    z = exp(-e*beta)/sum(exp(-e*beta)); % Partition function
    if temp == 0
        J0i(:,ii) = real(wavi(:,1,ii)'*Jz*wavi(:,1,ii));       
        J0j(:,ii) = real(wavj(:,1,ii)'*Jz*wavj(:,1,ii));
        Jsi(:,ii) = real(wav(:,1,ii)'*Jzhi*wav(:,1,ii));
        Jsj(:,ii) = real(wav(:,1,ii)'*Jzhj*wav(:,1,ii));
    else
        J0i(:,ii) = real(diag(wavi(:,:,ii)'*Jz*wavi(:,:,ii)))'*zi;
        J0j(:,ii) = real(diag(wavj(:,:,ii)'*Jz*wavj(:,:,ii)))'*zj;
        Jsi(:,ii) = real(diag(wav(:,:,ii)'*Jzhi*wav(:,:,ii)))'*z;
        Jsj(:,ii) = real(diag(wav(:,:,ii)'*Jzhj*wav(:,:,ii)))'*z;
    end
end

figure
plot(Hx,eni(2:end,:),'o')
hold on
plot(Hx,enj(2:end,:),'o')
plot(Hx,ens(2:4,:),'--k')
legend('sys1','sys2','s-s cpl')
xlabel('Magnetic field')
ylabel('Energy')

figure
plot(Hx,J0i(2,:),'-o')
hold on
plot(Hx,J0j(2,:),'-o')
plot(Hx,Jsi(2,:),'--k')
plot(Hx,Jsj(2,:),'--k')
legend('sys1','sys2','s-s cpl','s-s cpl')
xlabel('Magnetic field')
ylabel('<Jz>')
clearvars tempf intf scale