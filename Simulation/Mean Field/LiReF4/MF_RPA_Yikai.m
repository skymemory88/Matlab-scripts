function MF_linear_response(Temperatures,theta,phi,gama,mode)
% Current version assumes complete symmetrical equivalence among the four spin moments per unit cell
% TT: temperatures (can be an array). theta: angle between the external field and c-axis.
% phi: inplane angle in ab-plane between the magnetic field and a/b-axis. gamma: lifetime of hyperfine levels.
% scanMode: 1. Field plot with RPA. 2. wavevector plot with RPA

Options.RPA = true; % Apply random phase approximation (RPA) correction
Options.plotting = true; % Decide whether or not to plot the data at the end
Options.meV = true; % Energy unit choice (meV or GHz)
Options.saving = false; % Options to save the susceptibility tensors
Options.scanMode = mode; % 1. Field plot with RPA. 2. wavevector plot with RPA
clearvars -except Temperatures theta phi gama Options;

% Declare physical constants as global for consistency
global gLande_Ho ELEf NUCf dip_range muB J2meV mu0
muB = 9.274e-24; %[J/T]
mu0 = 4e-7*pi; % [H/m]
J2meV = 6.24151e+21; % Convert Joule to meV
dip_range = 100;
gLande_Ho = 1.25;
ELEf = gLande_Ho*muB*J2meV;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173

freq_total = (0:0.5:170);
% freq_total = (0:0.02:5);
for ii = 1:length(Temperatures)
    location = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\with Hz_I';
%     location = '/Volumes/GoogleDrive/My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I';
    filename = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg', Temperatures(ii), theta, phi),'.mat');
    file = fullfile(location,filename);
    load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" and "V"
    % which are eigenstates and eigenvalues calculated in the mean-field model
    % as a function of transverse field and temperature
    fields = vecnorm(fff,2,1);
    fprintf('Calculating for T = %.3f K.\n', Temperatures(ii));
    %         [fields, freq_total, rechi, imchi] = linear_response(eee,fff,ttt,vvv);
    if Options.RPA == true
        switch Options.scanMode %1. Field plot with RPA. 2. wavevector plot with RPA
            case 1
                qz = [0]';
%                 qx = [0.01 0.1 0.3 0.6 1]';
                qy = zeros(size(qz,1),1);
                qx = zeros(size(qz,1),1);
                qvec = [qx qy qz];
                eigenW = vvv;
                eigenE = eee;
            case 2
                qx = (1:0.01:2)';
                qy = zeros(size(qx,1),1);
                qz = zeros(size(qx,1),1);
                qvec = [qx qy qz];
                B0 = [4.3]; % choose the static magnetic field
                bidx = int16.empty(0,length(B0));
                for jj = 1:length(B0)
                    [~,bidx(jj)] = min(abs(vecnorm(fields,1,1)-B0(jj)));
                end
                fields = fff(:,bidx);
                eigenW = vvv(bidx,:,:);
                eigenE = eee(bidx,:);
            otherwise
                disp('Scan mode selection error!')
        end
        [fields, freq_total, chi0, ~] = linear_response(eigenE,fields,freq_total,ttt,eigenW,gama);
        [fields, freq_total, chiq.J] = RPA(qvec, fields, freq_total, ion, chi0.J);
%         [fields, freq_total, chiq.IJ] = RPA(qvec, fields, freq_total, ion, chi0.IJ);
%         [fields, freq_total, chiq.I] = RPA(qvec, fields, freq_total, ion, chi0.I);

%         chiq = ELEf^2*chiq.J + 2*ELEf*NUCf*chiq.IJ + NUCf^2*chiq.I;
%         chiq = ELEf^2*chiq.J + 2*ELEf*NUCf*chiq.IJ;
        chiq = ELEf^2*chiq.J;
    else
        [fields, freq_total, chi0, ~] = linear_response(eigenE,fields,freq_total,ttt,eigenW,gama);
    end
%   chi0 = ELEf^2*chi0.J + 2*ELEf*NUCf*chi0.IJ + NUCf^2*chi0.I;
    chi0 = ELEf^2*chi0.J;
    
    if Options.saving == true % Save the susceptibilities
        file_name1 = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e.mat', ttt, theta, phi,gama));
        savefile1 = fullfile(location,file_name1);
        if Options.RPA == true
            file_name2 = strcat('RPA_LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e.mat', ttt, theta, phi,gama));
            savefile2 = fullfile(location,file_name2);
            save_file(fields,freq_total,chi0,gama,savefile1); % Save the data without RPA corrections
            save_file(fields,freq_total,chiq,gama,savefile2); % save the data with RPA corrections
        else
            save_file(fields,freq_total,chi0,gama,savefile1); % Save the data without RPA corrections
        end
    end
    
    if Options.plotting == true % Plot the susceptibilities
        figs(fields,freq_total,chi0,gama,[0,0,0],Options,Temperatures(ii),'\chi_0^{zz}');
        if Options.RPA == true
            figs(fields,freq_total,chiq,gama,qvec,Options,Temperatures(ii),'\chi_{RPA}^{zz}');
            if Options.scanMode == 1
                figs(fields,freq_total,chiq-chi0,gama,qvec,1,Temperatures(ii),'\chi_{RPA}^{zz}-\chi_0^{zz}');
            end
        end
    end
end
end

function figs(fields,freq_total,chi,gama,qvec,Options,temperature,fig_tit)
meV = true;
f2E = 1/241.8;
if Options.meV == true
    freq_total = freq_total*f2E;
end
if Options.scanMode == 1
    for ii = 1:size(qvec,1)
        pos0 = [100 300 600 400]; % initial figure position
        pos_inc = [100 0 0 0];
%         % Color plot of the imaginary part of the susceptibility of x component
%         fig0 = figure;
%         set(fig0,'position',pos0);
%         hp0 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(1,1,:,:,ii))));
%         set(hp0, 'edgeColor','none')
%         caxis([0 5]);
%         colorbar
%         legend(['\gamma =' num2str(gama,'%.2e meV')]);
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         if meV == true
%             ylabel('Energy (meV)')
%         end
%         title(['Imaginary part of ', fig_tit, ' at ',sprintf('T = %.3f K.', temperature),sprintf('Q = [%1$.1f %2$.1f %3$.1f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])
% 
%         % Color plot of the imaginary part of the susceptibility of y component
%         fig1 = figure;
%         set(fig1,'position',pos0 + pos_inc);
%         hp1 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(2,2,:,:,ii)));
%         set(hp1, 'edgeColor','none')
%         caxis([0 5]);
%         colorbar
%         legend(['\gamma =' num2str(gama,'%.2e meV')]);
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         if meV == true
%             ylabel('Energy (meV)')
%         end
%         title(['Imaginary part of ', fig_tit, ' at ', sprintf('T = %.3f K.', temperature),sprintf('Q = [%1$.1f %2$.1f %3$.1f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])

        % Color plot the imaginary part of the susceptibilities of z component
        fig2 = figure;
        set(fig2,'position',pos0 + 2*pos_inc);
        hp2 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(3,3,:,:,ii))));
        set(hp2, 'edgeColor','none')
        caxis([0 5]);
        colorbar
        legend(['\gamma =' num2str(gama,'%.2e meV')]);
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        if meV == true
            ylabel('Energy (meV)')
        end
        title(['Imaginary part of ', fig_tit, ' at ',sprintf('T = %.3f K and ', temperature),sprintf('Q = [%1$.1f %2$.1f %3$.1f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])

        % Plot the real part of the susceptibility of the z component
        fig3 = figure;
        set(fig3,'position',pos0 + 3*pos_inc);
        hp3 = pcolor(fields(1,:),freq_total,squeeze(real(chi(3,3,:,:,ii))));
        set(hp3, 'edgeColor','none')
        caxis([0 5]);
        colorbar
        legend(['\gamma =' num2str(gama,'%.2e meV')]);
        xlabel('Magnetic field (T)')
        ylabel('Frequency (GHz)')
        if meV == true
            ylabel('Energy (meV)')
        end
        title(['Real part of ', fig_tit, ' at ',sprintf('T = %.3f K and ', temperature),sprintf('Q = [%1$.1f %2$.1f %3$.1f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])

    % % Plot the expectation value of <Jz+Iz>
    %     fig4 = figure;
    %     set(fig4,'position',pos0 + 4*pos_inc);
    %     hp4 = plot(fields(1,:),squeeze(Jz_exp),'-o');
    %     xlabel('Magnetic field (T)')
    %     ylabel('<J_z + I_z>')
    %     title({'Expectation value of J_z + I_z', 'in z direction'})
    end
elseif Options.scanMode == 2
    qvec = vecnorm(qvec,1,2);
    for ii = 1:size(fields,2)
        pos0 = [100 300 600 400]; % initial figure position
        pos_inc = [100 0 0 0];
        % Color plot of the imaginary part of the susceptibility of x component
%         fig0 = figure;
%         set(fig0,'position',pos0);
%         hp0 = pcolor(qvec,freq_total,(squeeze(imchi.x(:,ii,:))));
%         set(hp0, 'edgeColor','none')
%         set(gca, 'xdir', 'reverse' )
%         caxis([-10 10]);
%         colorbar
%         legend(['\gamma =' num2str(gama,'%.2e meV')]);
%         xlabel(sprintf('Q = [h, 0, 0]'))
%         ylabel('Frequency (GHz)')
%         if meV == true
%             ylabel('Energy (meV)')
%         end
%         title(['Imaginary part of ', fig_tit, ' at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])
% 
%         % Color plot of the imaginary part of the susceptibility of y component
%         fig1 = figure;
%         set(fig1,'position',pos0 + pos_inc);
%         hp1 = pcolor(qvec,freq_total,(squeeze(imchi.y(:,ii,:))));
%         set(hp1, 'edgeColor','none')
%         set(gca, 'xdir', 'reverse' )
%         caxis([0 5]);
%         colorbar
%         legend(['\gamma =' num2str(gama,'%.2e meV')]);
%         xlabel(sprintf('Q = [h, 0, 0]'))
%         ylabel('Frequency (GHz)')
%         if meV == true
%             ylabel('Energy (meV)')
%         end
%         title(['Imaginary part of ', fig_tit, ' at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])

        % Color plot the imaginary part of the susceptibilities of z component
        fig2 = figure;
        set(fig2,'position',pos0 + 2*pos_inc);
        hp2 = pcolor(qvec,freq_total,real(squeeze(chi(3,3,:,ii,:))));
        set(hp2, 'edgeColor','none')
        set(gca, 'xdir', 'reverse' )
        caxis([0 100]);
        colorbar
        legend(['\gamma =' num2str(gama,'%.2e meV')]);
        xlabel(sprintf('Q = [h, 0, 0]'))
        ylabel('Frequency (GHz)')
        if meV == true
            ylabel('Energy (meV)')
        end
        title(['Imaginary part of ', fig_tit, ' at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])
        
        % Plot the real part of the susceptibility of the z component
        fig3 = figure;
        set(fig3,'position',pos0 + 3*pos_inc);
        hp3 = pcolor(qvec,freq_total,imag(squeeze(chi(3,3,:,ii,:))));
        set(hp3, 'edgeColor','none')
        set(gca, 'xdir', 'reverse' )
        caxis([0 100]);
        colorbar
        legend(['\gamma =' num2str(gama,'%.2e meV')]);
        xlabel(sprintf('Q = [h, 0, 0]'))
        ylabel('Frequency (GHz)')
        if meV == true
            ylabel('Energy (meV)')
        end
        title(['Real part of ', fig_tit, ' at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])
    end
else
    disp('Unknow scan mode!')
end
end

function save_file(fields,freq_total,chi,gama,savefile)

% save(strcat('LiHoF4_x1y_x2y_',sprintf('%1$3.3fK_%2$uDeg', ttt, theta),'fields','freq_total','chi','gama'));
if size(chi,5) == 1
    x1z = squeeze(real(chi(3,3,:,:)));
    x2z = squeeze(imag(chi(3,3,:,:)));
    save(savefile,'fields','freq_total','x1z','x2z','gama','-v7.3');
else
    x1z = squeeze(real(chi(3,3,:,:,:)));
    x2z = squeeze(imag(chi(3,3,:,:,:)));
    save(savefile,'fields','freq_total','x1z','x2z','gama','-v7.3');
end
end

function [fields, freq_total, chi0, JIz_exp]=linear_response(eigenE,fields,freq_total,temperature,eigenW,gama)
global ELEf NUCf  
% Calculation of susceptibility tensor
chi0r_J = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0i_J = zeros(3,3,length(freq_total(1,:)),size(fields,2));
% chi0r_I = zeros(3,3,length(freq_total(1,:)),size(fields,2));
% chi0i_I = zeros(3,3,length(freq_total(1,:)),size(fields,2));
% chi0r_IJ = zeros(3,3,length(freq_total(1,:)),size(fields,2));
% chi0i_IJ = zeros(3,3,length(freq_total(1,:)),size(fields,2));

%Initiate J operators
J = 8; % Ho
I = 3.5; % Ho
Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
JhT.z=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*I+1)); % Expand Hilbert space
Jmh=kron(Jm,eye(2*I+1)); % Expand Hilbert space
JhT.x=(Jph+Jmh)/2;
JhT.y=(Jph-Jmh)/2i;

%Initiate I operators
Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
IhT.z=kron(eye(2*J+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*J+1),Ip); % Expand to match the Hilbert space
Imh=kron(eye(2*J+1),Im); % Expand to match the Hilbert space
IhT.x=(Iph+Imh)/2;
IhT.y=(Iph-Imh)/2i;
        
psdJhT.x = JhT.x+NUCf/ELEf*IhT.x;
psdJhT.y = JhT.y+NUCf/ELEf*IhT.y;
psdJhT.z = JhT.z+NUCf/ELEf*IhT.z;

% Single out <Jz+Iz> calculations
JIz_exp = double.empty(size(fields,2),size(JhT.z,1),0); % Expectation value of J-I pseudo-spin
for kk = 1:size(fields,2) % calculate susceptibility for all fields
    v = squeeze(eigenW(kk,:,:)); % Obtain the corresponding eigen vectors
    JzhT = JhT.z * ELEf;
    IzhT = IhT.z * NUCf;
%     JzhT = JhT.z;
%     IzhT = 0;
    ttz  = v'  * (JzhT+IzhT) * v;
    JIz_exp(kk,:,1) = real(diag(ttz));
end

% figure
% plot(fields, squeeze(JIz_exp(:,1:8,:)),'--')
% xlabel('Magnetic field (T)')
% ylabel('<J\times I>')
% ylabel('<J>')

for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (m);
    f2E = 1/241.8;  % GHz to meV
    omega = freq*f2E;   % define frequency sweep range (meV)
%     for k = 1:size(fields,2) % for debugging: calculate susceptibility for all fields
    parfor k = 1:size(fields,2) % calculate susceptibility for all fields
        v = squeeze(eigenW(k,:,:)); % Obtain the corresponding eigen vectors
        en = squeeze(eigenE(k,:)); % Obtain the corresponding eigen energies in meV
        if temperature ~= 0
            beta = 11.6/temperature; %[meV^-1]
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
        [ee,eep] = meshgrid(en,en);
        EE = eep-ee-omega;
        gamma = ones(size(EE))*gama;
        G = gamma ./ (EE.^2 + gamma.^2);
        G1 = EE ./ (EE.^2 + gamma.^2);

        [chi0r_J(:,:,m,k), chi0i_J(:,:,m,k)] = chi_Mx(psdJhT,psdJhT,v,NN,G,G1);
%         [chi0r_J(:,:,m,k), chi0i_J(:,:,m,k)] = chi_Mx(JhT,JhT,v,NN,G,G1);
%         [chi0r_IJ(:,:,m,k), chi0i_IJ(:,:,m,k)] = chi_Mx(IhT,JhT,v,NN,G,G1);
%         [chi0r_I(:,:,m,k), chi0i_I(:,:,m,k)] = chi_Mx(IhT,IhT,v,NN,G,G1);
    end
end

chi0.J = (chi0r_J + 1i*chi0i_J); % Electronic susceptibilities
% chi0.I = (chi0r_I + 1i*chi0i_I); % Nuclear susceptibilities
% chi0.IJ = (chi0r_IJ + 1i*chi0i_IJ); % Electronuclear Cross term

% chi0 = chi0_J + chi0_I + chi0_IJ;
end

function [xr,xi] = chi_Mx(opr1, opr2, wav, NN, G, G1)
% Calculation of matrix elements for the real and imaginary susceptibilities
% opr1: operator 1
% opr2: operator 2
% wav: eigenfunction of the hamiltonian
% NN: population differene between two eigenstates
% G: prefactor for real part of the susceptibility
% G1: prefactor for imaginary part of the susceptibility

ttx1  = wav'  * opr1.x * wav;
tty1  = wav'  * opr1.y * wav;
ttz1  = wav'  * opr1.z * wav;

ttx2  = wav'  * opr2.x * wav;
tty2  = wav'  * opr2.y * wav;
ttz2  = wav'  * opr2.z * wav;

% Calculate susceptibilities along a-axis
chi_tx  = (ttx1) .* (ttx2.') .* NN .* G;
chi_t1x = (ttx1) .* (ttx2.') .* NN .* G1;
xx = sum(sum(chi_tx));
xx1 = sum(sum(chi_t1x));

% Calculate susceptibilities along b-axis
chi_ty  = (tty1) .* (tty2.') .* NN .* G;
chi_t1y = (tty1) .* (tty2.') .* NN .* G1;
yy = sum(sum(chi_ty));
yy1 = sum(sum(chi_t1y));

% Calculate susceptibilities along c-axis
chi_tz  = (ttz1) .* (ttz2.') .* NN .* G;
chi_t1z = (ttz1) .* (ttz2.') .* NN .* G1;
zz = sum(sum(chi_tz));
zz1 = sum(sum(chi_t1z));

% Calculate susceptibilities along ab-axis
chi_txy  = (ttx1) .* (tty2.') .* NN .* G;
chi_t1xy = (ttx1) .* (tty2.') .* NN .* G1;
xy = sum(sum(chi_txy));
xy1 = sum(sum(chi_t1xy));

% Calculate susceptibilities along ac-axis
chi_txz  = (ttx1) .* (ttz2.') .* NN .* G;
chi_t1xz = (ttx1) .* (ttz2.') .* NN .* G1;
xz = sum(sum(chi_txz));
xz1 = sum(sum(chi_t1xz));

% Calculate susceptibilities along ba-axis
chi_tyx  = (tty1) .* (ttx2.') .* NN .* G;
chi_t1yx = (tty1) .* (ttx2.') .* NN .* G1;
yx = sum(sum(chi_tyx));
yx1 = sum(sum(chi_t1yx));

% Calculate susceptibilities along bc-axis
chi_tyz  = (tty1) .* (ttz2.') .* NN .* G;
chi_t1yz = (tty1) .* (ttz2.') .* NN .* G1;
yz = sum(sum(chi_tyz));
yz1 = sum(sum(chi_t1yz));

% Calculate susceptibilities along ca-axis
chi_tzx  = (ttz1) .* (ttx2.') .* NN .* G;
chi_t1zx = (ttz1) .* (ttx2.') .* NN .* G1;
zx = sum(sum(chi_tzx));
zx1 = sum(sum(chi_t1zx));

% Calculate susceptibilities along cb-axis
chi_tzy  = (ttz1) .* (tty2.') .* NN .* G;
chi_t1zy = (ttz1) .* (tty2.') .* NN .* G1;
zy = sum(sum(chi_tzy));
zy1 = sum(sum(chi_t1zy));

xr = [xx1 xy1 xz1
    yx1 yy1 yz1
    zx1 zy1 zz1];

xi = [xx xy xz
    yx yy yz
    zx zy zz];
end

function [fields, freq_total, chiq] = RPA(qvec, fields, freq_total, ion, chi0)

global muB mu0 J2meV gLande_Ho dip_range
N = 4; % Number of magnetic atoms in unit cell
a = [5.175 0 0;
     0 5.175 0;
     0 0 10.75]; % Lattice constant for LiHoF4
gfac = (muB)^2*(mu0/4/pi)*J2meV*10^30;

ex_range = 1; % Exchange interaction range (number of unit cell)
chiq = zeros(3,3,length(freq_total(1,:)),size(fields,2),size(qvec,1));
% chiqz = zeros(length(freq_total(1,:)),size(fields,2),size(qvec,1));
D = zeros(3,3,N,N,size(qvec,1));
for jj=1:size(qvec,1)
    D(:,:,:,:,jj) = gLande_Ho^2*(gfac*dipole_direct(qvec(jj,:),dip_range,a))...
        + exchange(qvec(jj,:),ion.ex(2),a,ex_range);
end

deno = zeros(3,3,size(freq_total,1),size(fields,2)); % RPA correction factor (denominator)
for ii = 1:size(fields,2)
    for kk = 1:length(freq_total(1,:))
        for nq = 1:size(qvec,1)
            MM = sum(sum(D(:,:,:,:,nq),4),3)/4; % average over the four sites in the unit cell and avoid double counting
%             chiqz(kk,ii,nq) = (1-MM(3,3)*chi0(3,3))\chi0(3,3);
            MM = chi0(:,:,kk,ii).*MM;
%             deno(:,:,kk,ii) = (eye(size(MM))- MM);
            deno(:,:,kk,ii) = (ones(size(MM))- MM);
            chiq(:,:,kk,ii,nq) = squeeze(deno(:,:,kk,ii)).\chi0(:,:,kk,ii);
%             chiq(:,:,kk,ii,nq) = squeeze(deno(:,:,kk,ii))\chi0(:,:,kk,ii);
        end
    end
end

% figure
% hp0 = pcolor(fields,freq_total,real(squeeze(chiqz)));
% set(hp0, 'edgeColor','none');
% colorbar
% xlabel('Magnetic field (T)')
% ylabel('Frequency (GHz)')
% title('Real part of \chi_{zz}^{RPA}')
% figure
% hp0 = pcolor(fields,freq_total,imag(squeeze(chiqz)));
% set(hp0, 'edgeColor','none');
% colorbar
% xlabel('Magnetic field (T)')
% ylabel('Frequency (GHz)')
% title('Imaginary part of \chi_{zz}^{RPA}')

% for ii = 1:3
%     figure
%     hp0 = pcolor(fields,freq_total,squeeze((Mx(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
% 
%     figure
%     hp0 = pcolor(fields,freq_total,squeeze(abs(deno(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
% end
end