function [ion,evolution,E,V,h_mf1] = remf(ion,hvec,t,withdemagn,alpha)
global strategies muB mu0 J2meV rundipole;
persistent dipole;

for ii = 1:size(ion.name,1)
    if ion.prop(ii)~=0
        elem = ii;
        if ion.prop(ii)==1
            if ion.hyp(ii) > 0
                H_dims = (2*ion.I(elem)+1)*(2*ion.J(elem)+1); % Calculate the dimension of the Hilbert space
            elseif ion.hyp(ii) == 0
                H_dims = 2*ion.J(elem)+1; % Calculate the dimension of the Hilbert space
            end
        end
    end
end

dip_range = 100; %The range over which the dipolar interaction is included
nn = 1; % range of exchange interaction in unit of unit cell dimension
if(isempty(dipole))
    dipole = dipole_direct([0 0 0],dip_range,ion.abc{elem});
    for ii = 1:size(ion.name,1)
        ion.cf(:,:,ii)={cf(ion.J(ii),ion.B(ii,:),ion.cfRot(ii))};
        ion.exch(:,:,ii)={exchange([0,0,0],ion.ex(ii),ion.abc{elem},nn)};
    end
end
if rundipole == true
    dipole = dipole_direct([0 0 0],dip_range,ion.abc{elem});
    for ii=1:size(ion.name,1)
        ion.cf(:,:,ii)={cf(ion.J(ii),ion.B(ii,:),ion.cfRot(ii))};
        ion.exch(:,:,ii)={exchange([0,0,0],ion.ex(ii),ion.abc{elem},nn)};
    end
    rundipole = false;
end

alt = [1  1  1
      -1  1  1
      -1 -1  1
       1 -1  1]; % Matrix for alternating summation over unit cell
   
% alt = [1  1  1
%       -1  1  1
%       -1 -1  1
%        1 -1  1]; % Original code (2021.06.20) --Yikai

momente_mean = 0;
ion.mom_hyp = ion.mom; % For virtual crystal approximation for the isotops
for ii=1:size(ion.name,1)
%     momente_mean = momente_mean + ion.prop(i)*ion.mom(:,:,i); % Original code (Yikai) --2021.06.20
    momente_mean = momente_mean + ion.prop(ii)*((1-ion.hyp(ii))*ion.mom(:,:,ii) + ion.hyp(ii)*ion.mom_hyp(:,:,ii));
end

eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:) = 1;
demagn = zeros(3,3,4,4);
demagn_t = ellipsoid_demagn(alpha);
% demagn_t = [7.55 0 0; 0 1.68 0; 0 0 -0.32]; % Thesis "Ultra-low temperature dilatometry" by John L. Dunn
demagn(1,1,:,:) = demagn_t(1,1);
demagn(2,2,:,:) = demagn_t(2,2);
demagn(3,3,:,:) = demagn_t(3,3);

Vc = sum(ion.abc{elem}(1,:).*cross(ion.abc{elem}(2,:),ion.abc{elem}(3,:)))*10^-30; % Volume of unit cell [A^-3]
dmag_fac = 4/Vc * mu0/4/pi * muB^2 * J2meV; % demag_fac = N/V * mu_0/4pi * (mu_b)^2
% Lorenz = dmag_fac * 4*pi/3; % Lorenz_term = N/V * mu_0/4pi * (mu_b)^2 * 4pi/3
% Lorenz = 3.124032/1000; % Original code (Ho)
% Lorenz = 3.15452/1000; % Original code (Er)
gfac = mu0*muB^2/4/pi*J2meV*10^30; % [meV*Ang^-3] (dipole_direct gives [Ang^3])

ion.d = zeros([3,3,4,4,size(ion.name,1)]);
ion.d_ex = zeros([3,3,4,4,size(ion.name,1)]);
for ii=1:size(ion.name,1)
%     ion.d(:,:,:,:,i) = ion.gLande(ii)*(gfac*dipole + eins*Lorenz/4 - withdemagn*dmag_fac*pi*demagn);
    ion.d(:,:,:,:,ii) = ion.gLande(ii)*(gfac*dipole + dmag_fac*pi*(eins/3 - withdemagn*demagn));
    ion.d_ex(:,:,:,:,ii) = ion.exch{:,:,ii};
    evolution.(ion.name{ii})(1,:,:) = ion.mom(:,:,ii);
    evolution.(ion.name_hyp{ii})(1,:,:) = ion.mom_hyp(:,:,ii);
end

% set convergence criteria. Convergence is accepted when either:
% (Difference < ConvMin) or (Difference < ConvMax and Niter>NiterMin) or (Niter>NiterMax)
ConvMax = 1e-5;
ConvMin = 1e-6;
NiterMin = 500;
NiterMax = 20000;

ion.mom_old = zeros([4,3,size(ion.name,1)]);
ion.mom_old_hyp = zeros([4,3,size(ion.name,1)]);

for iterations=1:NiterMax
    for ii=1:size(ion.name,1)
        ion.mom_old(:,:,ii) = ion.mom(:,:,ii);
        ion.mom_old_hyp(:,:,ii) = ion.mom_hyp(:,:,ii);
    end
    momente_old = momente_mean;
%     E=zeros(1,length(energies(1,:)));
    energies = zeros(4, H_dims); % H_dims = (2*I+1) x (2*J+1)
%     energies = zeros(4, 2); % truncate the energy levels to Ising basis
    for ionn=1:4 % each Ho3+ ion in the unit cell
        % calculate MF magnetic moments
        h_dipol = zeros([size(ion.name,1),3]);
        h_dipol_hyp = zeros([size(ion.name,1),3]);
        h_ex = zeros([size(ion.name,1),3]);
        h_ex_hyp = zeros([size(ion.name,1),3]);
        for ionm=1:4 
            for ii=1:size(ion.name,1)
                h_dipol(ii,:) = h_dipol(ii,:) + ion.mom_old(ionm,:,ii) * diag(ion.renorm(ii,:)) * ion.d(:,:,ionm,ionn,ii)';
                h_ex(ii,:) = h_ex(ii,:) + ion.mom_old(ionm,:,ii) * diag(ion.renorm(ii,:)) * ion.d_ex(:,:,ionm,ionn,ii)';
                h_dipol_hyp(ii,:) = h_dipol_hyp(ii,:) + ion.mom_old_hyp(ionm,:,ii) * diag(ion.renorm(ii,:)) * ion.d(:,:,ionm,ionn,ii)';
                h_ex_hyp(ii,:) = h_ex_hyp(ii,:) + ion.mom_old_hyp(ionm,:,ii) * diag(ion.renorm(ii,:)) * ion.d_ex(:,:,ionm,ionn,ii)';
            end
        end

        % Virtual meanfield
        h_mf = zeros([size(ion.name,1),3]);
        for ii=1:size(ion.name,1)
            h_mf(ii,:) = ion.prop(ii) * ((1-ion.hyp(ii)) * ( ion.gLande(ii)*h_dipol(ii,:) + h_ex(ii,:) )...
                + ion.hyp(ii) * (ion.gLande(ii) * h_dipol_hyp(ii,:) + h_ex_hyp(ii,:)));
            for j = 1:size(ion.name,1)-1 % this may be problematic for doped samples --Yikai
                k = (ii+j>size(ion.name,1)) * size(ion.name,1); % cyclic (perioric) condition
                h_mf(ii,:) = h_mf(ii,:) + ion.prop(ii+j-k) * ((1-ion.hyp(ii+j-k)) * ion.gLande(ii) * h_dipol(ii+j-k,:) +...
                    ion.hyp(ii+j-k) * ion.gLande(ii) * h_dipol_hyp(ii+j-k,:)); % other kinds of ions
            end
        end

        %calculate moments of ions in a meanfield (ie diagonnalize the hamiltonian)
        for ii = 1:size(ion.name,1)
            if ion.prop(ii) > 0 % if the ion exists
                if ion.hyp(ii) == 0 % if there is no hyperfine interaction
                    [jx,jy,jz,energies(ionn,:),v] = MF_moments(hvec, h_mf(ii,:), t, ion.J(ii), ion.gLande(ii), ion.cf{:,:,ii}, ion.h4(ii));
                    ion.mom(ionn,:,ii) = update_moments([jx,jy,jz], evolution.(ion.name{ii}), ionn,iterations);
                else % if there is finite hyperfine interaction
                    [jx,jy,jz,energies(ionn,:),v] = MF_moments_hyper(hvec, h_mf(ii,:), t, ion.J(ii), ion.gLande(ii), ion.nLande(ii),...
                        ion.cf{:,:,ii}, ion.A(ii), ion.I(ii), ion.h4(ii));
                    ion.mom_hyp(ionn,:,ii) = update_moments([jx,jy,jz], evolution.(ion.name_hyp{ii}), ionn, iterations);
                end
                if ion.hyp(ii) ~= 1 % if there is only partial hyperfine interaction
                    [jx,jy,jz,~,~] = MF_moments(hvec, h_mf(ii,:),t,ion.J(ii),ion.gLande(ii),ion.cf{:,:,ii},ion.h4(ii));
                    ion.mom(ionn,:,ii) = update_moments([jx,jy,jz],evolution.(ion.name{ii}),ionn,iterations);
                end
             
            end
        end
%      E(ionn)=mean(energies(ionn,:)); % Why does it take an arithmetic average of all the energy levels? --Yikai (12.02.2020)
        
        if strategies.symmetry % HMR: copy configuration to symmetry equivalent sites --Yikai (10.02.2020)
            energies = repmat(energies(1,:),4,1);
            ion.mom_hyp = repmat(ion.mom_hyp(1,:,:),4,1,1);
            ion.mom = repmat(ion.mom(1,:,:),4,1,1);
            break % exits the for ionn=1:4 loop
        end 
    end

    momente_mean = 0;
    for ii=1:size(ion.name,1)
        momente_mean = momente_mean + ion.prop(ii)*((1-ion.hyp(ii))*ion.mom(:,:,ii) + ion.hyp(ii)*ion.mom_hyp(:,:,ii));
        evolution.(ion.name{ii})(size(evolution.(ion.name{ii}),1)+1,:,:) = ion.mom(:,:,ii);
        evolution.(ion.name_hyp{ii})(size(evolution.(ion.name_hyp{ii}),1)+1,:,:) = ion.mom_hyp(:,:,ii);
    end

    % Iterate untill convergence criterion reached
    % Added energy convergence check --Yikai (2020-05-07)
    JDiff = sum(sum(abs(momente_old-momente_mean)));
    EDiff = [energies(4,:)-energies(3,:);energies(3,:)-energies(2,:);energies(2,:)-energies(1,:);energies(1,:)-energies(4,:)];
    if (~any(abs(EDiff(:))>=ConvMin) && JDiff<ConvMin) || (JDiff<ConvMax && iterations>NiterMin) 
%         E = energies(1,:);
        E = mean(energies,1);
        V = v;
        h_mf1 = h_mf;
        if elem == 1 % Use alternating spin sum for Er
            disp(num2str([t,hvec,iterations,mean(alt.*ion.mom_hyp(:,:,elem))],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f)'))
        else
            disp(num2str([t,hvec,iterations,mean(momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <J>=(%3.3f,%3.3f,%3.3f)'))
        end
        evolution.iterations = iterations;
        return
    end
end
energies = squeeze(energies);
E = energies(1,:);
V = v;
h_mf1 = h_mf;
if elem == 1 % Use alternating spin sum for Er
    disp(num2str([t,hvec,iterations,mean(alt.*ion.mom_hyp(:,:,elem))],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f), fail to converge!'))
else
    disp(num2str([t,hvec,iterations,mean(momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <J>=(%3.3f,%3.3f,%3.3f), fail to converge!'))
end
evolution.iterations = iterations;
return

function J = update_moments(Jnew,evolution,ionn,iterations)
global strategies;
if strategies.accelerator~=0 && iterations>10 && mod(iterations,strategies.expfit_period)~=0 % use accellerated update to shorten convergence
    Jold = squeeze(evolution(end-1,ionn,:))';
    Jnow = squeeze(evolution(end,ionn,:))';
    dnow = Jnow-Jold;
    J = Jnew + strategies.accelerator*dnow;
elseif strategies.expfit && mod(iterations,strategies.expfit_period)==0
    ndif = strategies.expfit_deltaN;
    J = Jnew;
%     Jold=squeeze(evolution(end-2*ndif+1-4,ionn,:))'; Why the additional "-4"? --Yikai (12.02.2020)
%     Jnow=squeeze(evolution(end-1*ndif+1-4,ionn,:))';
%     Jnew=squeeze(evolution(end-0*ndif+1-4,ionn,:))';
    Jold = squeeze(evolution(end-2*ndif+1,ionn,:))';
    Jnow = squeeze(evolution(end-1*ndif+1,ionn,:))';
    Jnew = squeeze(evolution(end-0*ndif,ionn,:))';
%    enNr=min(max(0,(Jnew-Jnow)./(Jnow-Jold)),0.998);
%    if isempty(find(enNr>=1,1))
    enNr = (Jnew-Jnow)./(Jnow-Jold);
    for nxyz=1:3 % for each of the three (dimensions) components of the Js --Yikai (10.02.2020)
       if enNr(nxyz)>0.01 && enNr(nxyz)<0.998
          feNNr = (Jnew-Jnow)./(enNr-1);
          J(nxyz) = J(nxyz)-feNNr(nxyz);
       end
    end
%    if isempty(find(enNr>=1,1))
%        feNNr = (Jnew-Jnow)./(enNr-1);
%        J = Jnow-feNNr;
%    else
%        J = Jnew;
%    end
else
    J = strategies.damping*squeeze(evolution(end,ionn,:))'+(1-strategies.damping)*Jnew;
end
return

function [jx,jy,jz,energies,wv] = MF_moments(hvec,h_int,t,J,gLande,Hcf,h4)
global muB J2meV
ELEf = gLande*muB*J2meV;     % Lande factor * Bohr magneton (meV T^-1)

% Initiate J operators
Jz = diag(J:-1:-J);
Jp = diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm = Jp';
Jx = (Jp+Jm)/2;
Jy = (Jp-Jm)/2i;

% Calculate Hamiltonian
Hzeeman = -ELEf*(hvec(1)*Jx + hvec(2)*Jy + hvec(3)*Jz); % Electronic Zeeman interaction [meV]
Hint = -(h_int(1)*Jx + h_int(2)*Jy + h_int(3)*Jz); % dipole-dipole interaction
% Hint = -h_int(3)*Jz; % include only diagonal interaction
O44 = (Jp^4+Jm^4)/2; % off-diagonal CEF term for in-plane anisotropy
H_h4 = h4*O44*h_int(1)^2; % order by disorder anisotropy

Ham = Hcf + Hzeeman + Hint + H_h4;
% Ham = Hcf; % Crystal electrostatic field terms only
% Ham = Hcf + Hzeeman + H_h4; % without dipole interaction
% Ham = Hzeeman + HvM + H_h4; % without crystal field term


% Diagonalize
[wv,ee] = eig(Ham); % [meV]
ee = real(diag(ee));
[ee,n] = sort(ee);
energies = ee;
ee = ee-min(ee); % normalize the eigenenergy to the ground state 

wv = wv(:,n); % reorder the eigenfunctions according to eigenstates
beta = 11.6/t; % [kB*T]^-1 in [meV]

% e = e(1:2); % truncate the vasis to form an Ising basis
% v = v(:,1:2); % truncate the basis to form an Ising basis

% Calculate Matrixelements and Moments
if t == 0 % At zero temperature, use only lowest eigenvalue.
    jx = real(wv(:,1)'*Jx*wv(:,1));
    jy = real(wv(:,1)'*Jy*wv(:,1));
    jz = real(wv(:,1)'*Jz*wv(:,1));
else % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    z = exp(-ee*beta)/sum(exp(-ee*beta)); % Density matrix element
    jx = real(diag(wv'*Jx*wv))'*z;
    jy = real(diag(wv'*Jy*wv))'*z;
    jz = real(diag(wv'*Jz*wv))'*z;
%     energies=sum(E)*z; % This assignment is very puzzling --Yikai
end

return

function [jx,jy,jz,energies,wv] = MF_moments_hyper(hvec,h_int,t,J,gLande,nLande,Hcf,A,I,h4)
% With hyperfine coupling
global muB muN J2meV Options
ELEf = gLande*muB*J2meV;     % Lande factor * Bohr magneton (meV/T)
NUCf = nLande * muN; % gN = mu/mu_N/<I> = 4.173/(7/2)
% NUCf = 1.668 * muN;  % http://easyspin.org/documentation/isotopetable.html

%Initiate J operators
Jz = diag(J:-1:-J);
Jzh = kron(Jz,eye(2*I+1));
Jp = diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm = Jp';
Jph = kron(Jp,eye(2*I+1));
Jmh = kron(Jm,eye(2*I+1));
Jxh = (Jph+Jmh)/2;
Jyh = (Jph-Jmh)/2i;

% Initiate I operators
Iz = diag(I:-1:-I);
Izh = kron(eye(2*J+1),Iz);
Ip = diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im = Ip';
Iph = kron(eye(2*J+1),Ip);
Imh = kron(eye(2*J+1),Im);
Ixh = (Iph+Imh)/2;
Iyh = (Iph-Imh)/2i;

Hcfh = kron(Hcf,eye(2*I+1)); % Expand the crystal field space to include nuclear moments' degrees of freedom

%Calculate Hamiltonian
O44 = (Jph^4+Jmh^4)/2; % off-diagonal CEF term for in-plane anisotropy
H_h4 = h4*O44*h_int(1)^2; % order by disorder anisotropy
Hint =  -(h_int(1)*Jxh + h_int(2)*Jyh + h_int(3)*Jzh); % virtual mean field
% Hint =  -h_int(3)*Jzh; % include only dipole interaction along z-axis
Hzeeman = -ELEf*(hvec(1)*Jxh + hvec(2)*Jyh + hvec(3)*Jzh); % Electron Zeeman interaction
Hyper = A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh); % hyperfine interaction

if Options.nZee == true 
    HzeemanI = -NUCf*(hvec(1)*Ixh + hvec(2)*Iyh + hvec(3)*Izh); % Nuclear Zeeman interaction
    Ham = Hcfh + H_h4 + Hint + Hzeeman + HzeemanI + Hyper;
    % Ham = Hcfh + H_h4 + Hzeeman + HzeemanI + A*(Izh*Jzh); % with Ising hyperfine interaction
else
%     Ham = Hzeeman + Hyper + HvM; % barebone hamiltonian
    Ham = Hcfh + H_h4 + Hint + Hzeeman + Hyper;
    % Ham = Hcfh + H_h4 + Hzeeman + A*(Izh*Jzh); % with Ising hyperfine interaction
end

% introduce order by disorder anisotropy
O44c = (Jph^4+Jmh^4)/2;
Ham = Ham + h4*O44c*h_int(1)^2;

% Diagonalize (Original)
[wv,ee] = eig(Ham);
ee = real(diag(ee)); % Take only the real part of the eigen-energy to form a diaganol matrix
[ee, n] = sort(ee); % sort the energy from lowest to the highest

energies = ee; % save the unnormalized eigenenergies
ee = ee-min(ee); % Normalize the energy amplitude to the lowest eigen-energy

wv = wv(:,n); % sort the eigen-vectors in its basis accordingly
beta = 11.6/t; % [kB*T]^-1 in [meV]

% Calculate Matrix elements and Moments
if t==0 
    % energies=e(1,1); % At zero temperature, use only lowest eigenvalue.
    jx = real(wv(:,1)'*Jxh*wv(:,1));
    jy = real(wv(:,1)'*Jyh*wv(:,1));
    jz = real(wv(:,1)'*Jzh*wv(:,1));
else % Boltzman factor (with t in Kelvin)
    z = exp(-ee*beta)/sum(exp(-ee*beta)); % Partition function
    jx = real(diag(wv'*Jxh*wv))'*z; % Calculate the expectation values
    jy = real(diag(wv'*Jyh*wv))'*z;
    jz = real(diag(wv'*Jzh*wv))'*z;
%     energies = sum(E)*z; % this is a puzzling assignment (Yikai)
end

return