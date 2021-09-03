function [ion,evolution,E,V,h_mf1] = remf(ion,hvec,t,withdemagn,alpha)
global strategies muB mu0 J2meV rundipole;
persistent dipole;

for i=1:size(ion.name,1)
    if ion.prop(i)~=0
        elem = i;
        if ion.prop(i)==1
            if ion.hyp(i) > 0
                H_dims = (2*ion.I(elem)+1)*(2*ion.J(elem)+1); % Claculate the dimension of the Hilbert space
            elseif ion.hyp(i) == 0
                H_dims = 2*ion.J(elem)+1; % Claculate the dimension of the Hilbert space
            end
        end
    end
end

dip_range = 100; %The range over which the dipolar interaction is included
nn = 1; % range of exchange interaction in unit of unit cell dimension
if(isempty(dipole))
    dipole=dipole_direct([0 0 0],dip_range,ion.abc{elem});
    for i=1:size(ion.name,1)
        ion.cf(:,:,i)={cf(ion.J(i),ion.B(i,:))};
        ion.exch(:,:,i)={exchange([0,0,0],ion.ex(i),ion.abc{elem},nn)};
    end
end
if rundipole == true
    dipole=dipole_direct([0 0 0],dip_range,ion.abc{elem});
    for i=1:size(ion.name,1)
        ion.cf(:,:,i)={cf(ion.J(i),ion.B(i,:))};
        ion.exch(:,:,i)={exchange([0,0,0],ion.ex(i),ion.abc{elem},nn)};
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
for i=1:size(ion.name,1)
%     momente_mean = momente_mean + ion.prop(i)*ion.mom(:,:,i); % Original code (Yikai) --2021.06.20
    momente_mean = momente_mean + ion.prop(i)*((1-ion.hyp(i))*ion.mom(:,:,i) + ion.hyp(i)*ion.mom_hyp(:,:,i));
end

%gLande
for i=1:size(ion.name,1)
    ion.gLande(i)=gLande(ion.L(i),ion.S(i));
end

eins = zeros(3,3,4,4);
eins(1,1,:,:) = 1; eins(2,2,:,:) = 1; eins(3,3,:,:)=1;
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

ion.d=zeros([3,3,4,4,size(ion.name,1)]);
ion.d_ex=zeros([3,3,4,4,size(ion.name,1)]);
for i=1:size(ion.name,1)
%     ion.d(:,:,:,:,i) = ion.gLande(i)*(gfac*dipole + eins*Lorenz/4 - withdemagn*dmag_fac*pi*demagn);
    ion.d(:,:,:,:,i) = ion.gLande(i)*(gfac*dipole + dmag_fac*pi*(eins/3 - withdemagn*demagn));
    ion.d_ex(:,:,:,:,i) = ion.exch{:,:,i};
    evolution.(ion.name{i})(1,:,:) = ion.mom(:,:,i);
    evolution.(ion.name_hyp{i})(1,:,:) = ion.mom_hyp(:,:,i);
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
    for i=1:size(ion.name,1)
        ion.mom_old(:,:,i) = ion.mom(:,:,i);
        ion.mom_old_hyp(:,:,i) = ion.mom_hyp(:,:,i);
    end
    momente_old = momente_mean;
%     E=zeros(1,length(energies(1,:)));
    energies = zeros(4, H_dims); % H_dims = (2*I+1) x (2*J+1)
%     energies = zeros(4,2); % truncate the energy levels to Ising basis
    for ionn=1:4 % each Ho3+ ion in the unit cell
        %calculate meanfield
        h_dipol = zeros([size(ion.name,1),3]);
        h_dipol_hyp = zeros([size(ion.name,1),3]);
        h_ex = zeros([size(ion.name,1),3]);
        h_ex_hyp = zeros([size(ion.name,1),3]);
        for ionm=1:4 
            for i=1:size(ion.name,1)
                h_dipol(i,:) = h_dipol(i,:) + ion.mom_old(ionm,:,i)*diag(ion.renorm(i,:))*ion.d(:,:,ionm,ionn,i)';
                h_ex(i,:) = h_ex(i,:) + ion.mom_old(ionm,:,i)*diag(ion.renorm(i,:))*ion.d_ex(:,:,ionm,ionn,i)';
                h_dipol_hyp(i,:) = h_dipol_hyp(i,:) + ion.mom_old_hyp(ionm,:,i)*diag(ion.renorm(i,:))*ion.d(:,:,ionm,ionn,i)';
                h_ex_hyp(i,:) = h_ex_hyp(i,:) + ion.mom_old_hyp(ionm,:,i)*diag(ion.renorm(i,:))*ion.d_ex(:,:,ionm,ionn,i)';
            end
        end

        % Virtual meanfield
        h_mf=zeros([size(ion.name,1),3]);
        for i=1:size(ion.name,1)
            h_mf(i,:) = ion.prop(i)*((1-ion.hyp(i))*(ion.gLande(i)*h_dipol(i,:)+h_ex(i,:)) + ion.hyp(i)*(ion.gLande(i)*h_dipol_hyp(i,:)+h_ex_hyp(i,:)));
            for j = 1:size(ion.name,1)-1 % this may be problematic for doped samples --Yikai
                k = (i+j>size(ion.name,1))*size(ion.name,1); % cyclic (periorid) condition
                h_mf(i,:) = h_mf(i,:) + ion.prop(i+j-k)*((1-ion.hyp(i+j-k))*ion.gLande(i)*h_dipol(i+j-k,:) +...
                    ion.hyp(i+j-k)*ion.gLande(i)*h_dipol_hyp(i+j-k,:));% other kinds of ions
            end
        end

        %calculate moments of ions in a meanfield (ie diagonnalize the hamiltonian)
        for i=1:size(ion.name,1)
            if ion.prop(i) > 0 % if the ion exists
                if ion.hyp(i) > 0 % if there is finite hyperfine interaction
                 [jx,jy,jz,energies(ionn,:),v] = MF_moments_hyper(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i},ion.A(i),ion.I(i),ion.h4);
                 ion.mom_hyp(ionn,:,i) = update_moments([jx,jy,jz],evolution.(ion.name_hyp{i}),ionn,iterations);
                end
                if ion.hyp(i) ~= 1 % if there is partial or no hyperfine interaction
                 [jx,jy,jz,~,~] = MF_moments(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i},ion.h4);
%                  [jx,jy,jz,energies(ionn,:),v] = MF_moments(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i},ion.h4);
                 ion.mom(ionn,:,i) = update_moments([jx,jy,jz],evolution.(ion.name{i}),ionn,iterations);
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
    for i=1:size(ion.name,1)
        momente_mean = momente_mean + ion.prop(i)*((1-ion.hyp(i))*ion.mom(:,:,i) + ion.hyp(i)*ion.mom_hyp(:,:,i));
        evolution.(ion.name{i})(size(evolution.(ion.name{i}),1)+1,:,:) = ion.mom(:,:,i);
        evolution.(ion.name_hyp{i})(size(evolution.(ion.name_hyp{i}),1)+1,:,:) = ion.mom_hyp(:,:,i);
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
    J = Jnew+strategies.accelerator*dnow;
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
%        feNNr=(Jnew-Jnow)./(enNr-1);
%        J=Jnow-feNNr;
%    else
%        J=Jnew;
%    end
else
    J = strategies.damping*squeeze(evolution(end,ionn,:))'+(1-strategies.damping)*Jnew;
end
return

function [jx,jy,jz,energies,v] = MF_moments(hvec,h_dipol,t,J,gLande,Hcf,h4)
global muB J2meV
ELEf = gLande*muB*J2meV;     % Lande factor * Bohr magneton (meV T^-1)

%Initiate J operators
Jz = diag(J:-1:-J);
Jp = diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm = Jp';
Jx = (Jp+Jm)/2;
Jy = (Jp-Jm)/2i;

%Calculate Hamiltonian
Hzeeman = -ELEf*(hvec(1)*Jx + hvec(2)*Jy + hvec(3)*Jz); %in meV (Yikai)
HvM = -(h_dipol(1)*Jx+h_dipol(2)*Jy+h_dipol(3)*Jz);
Ham = Hcf + Hzeeman + HvM;

% introduce order by disorder anisotropy
O44 = (Jp^4+Jm^4)/2;
Ham = Ham + h4*O44*h_dipol(1)^2;

%Diagonalize
[v,e] = eig(Ham); %in meV
E = e;
e = real(diag(e));
e = e-min(e);
[e,n] = sort(e);
v = v(:,n);
% e = e(1:2); % truncate the vasis to form an Ising basis
% v = v(:,1:2); % truncate the basis to form an Ising basis
beta = 11.6/t; % [kB*T]^-1 in [meV]

%Calculate Matrixelements and Moments
if t==0 % At zero temperature, use only lowest eigenvalue.
    energies = E(1,1);
    jx = real(v(:,1)'*Jx*v(:,1));
    jy = real(v(:,1)'*Jy*v(:,1));
    jz = real(v(:,1)'*Jz*v(:,1));
else % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e = e-min(e);
    z = exp(-e*beta)/sum(exp(-e*beta)); % Density matrix element
    jx = real(diag(v'*Jx*v))'*z;
    jy = real(diag(v'*Jy*v))'*z;
    jz = real(diag(v'*Jz*v))'*z;
%     energies=sum(E)*z; % This assignment is very puzzling --Yikai
    energies = e;
end

return

function [jx,jy,jz,energies,v] = MF_moments_hyper(hvec,h_dipol,t,J,gLande,Hcf,A,I,h4)
% With hyperfine coupling
global muB muN J2meV
ELEf = gLande*muB*J2meV;     % Lande factor * Bohr magneton (meV/T)
NUCf = 1.192 * muN; % gN = mu/mu_N/<I> = 4.173/(7/2)
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
% Expand the crystal field tensor to include nuclear moments' degrees of freedom
Hcfh = kron(Hcf,eye(2*I+1));
%Initiate I operators
Iz = diag(I:-1:-I);
Izh = kron(eye(2*J+1),Iz);
Ip = diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im = Ip';
Iph = kron(eye(2*J+1),Ip);
Imh = kron(eye(2*J+1),Im);
Ixh = (Iph+Imh)/2;
Iyh = (Iph-Imh)/2i;

%Calculate Hamiltonian
Hzeeman = -ELEf*(hvec(1)*Jxh + hvec(2)*Jyh + hvec(3)*Jzh); % Electron Zeeman interaction
HzeemanI = -NUCf*(hvec(1)*Ixh + hvec(2)*Iyh + hvec(3)*Izh); % Nuclear Zeeman interaction
Hyper = A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh); % hyperfine interaction
HvM =  -(h_dipol(1)*Jxh + h_dipol(2)*Jyh + h_dipol(3)*Jzh); % virtual mean field
Ham = Hcfh + Hzeeman + HzeemanI  + Hyper + HvM;
% Ham = Hcfh + Hzeeman + Hyper + HvM;
% Ham = Hcfh + Hzeeman + HzeemanI  + Hyper + A*(Izh*Jzh);

% introduce order by disorder anisotropy
O44c = (Jph^4+Jmh^4)/2;
Ham = Ham + h4*O44c*h_dipol(1)^2;

% Diagonalize (Original)
[v,e] = eig(Ham);
e = real(diag(e)); % Take only the real part of the eigen-energy to form a diaganol matrix
[e,n] = sort(e); % sort the energy from lowest to the highest
v = v(:,n); % sort the eigen-vectors in its basis accordingly
e = e-min(e); % Normalize the energy amplitude to the lowest eigen-energy
beta = 11.6/t; % [kB*T]^-1 in [meV]

%Calculate Matrix elements and Moments
if t==0 
    energies = e;
    % energies=e(1,1); % At zero temperature, use only lowest eigenvalue.
    jx = real(v(:,1)'*Jxh*v(:,1));
    jy = real(v(:,1)'*Jyh*v(:,1));
    jz = real(v(:,1)'*Jzh*v(:,1));
else % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    z = exp(-e*beta)/sum(exp(-e*beta)); % Partition function
    jx = real(diag(v'*Jxh*v))'*z; % Calculate the expectation values
    jy = real(diag(v'*Jyh*v))'*z;
    jz = real(diag(v'*Jzh*v))'*z;
%     energies = sum(E)*z; % this is a puzzling assignment (Yikai)
    energies = e;
end

return