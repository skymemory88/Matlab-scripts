function  MF_RPA_Yikai(Temperatures,phi,gama,mode,saveopt)
% Current version assumes complete symmetrical equivalence among the four spin moments per unit cell
% Temperatures (can be an array).
% phi: inplane angle in ab-plane between the magnetic field and a/b-axis [degree]
% gama: linewidth of hyperfine levels [meV].
% scanMode: 'field': Magnetic Field scan. 'qvec': wavevector scan with RPA
% saveopt: true or false to save the susceptibility data

Options.RPA = true; % Apply random phase approximation (RPA) correction
Options.plotting = false; % Decide whether or not to plot the data at the end
Options.meV = false; % Energy unit choice (meV or GHz)
Options.saving = saveopt; % Options to save the susceptibility tensors
Options.scanMode = mode; % 1. Field plot with RPA. 2. wavevector plot with RPA
clearvars -except Temperatures theta phi gama Options;

% Declare physical constants as global for consistency
hbar = 1.055E-34; % Reduced Planck constant [J.s]
global gLande_Ho ELEf NUCf dip_range ex_range muB J2meV mu0 f2E
muB = 9.274e-24; %[J/T]
mu0 = 4e-7*pi; % [H/m]
J2meV = 6.24151e+21; % Joule to meV
f2E = hbar*2*pi*10^9*J2meV;% GHz to meV
dip_range = 100; % dipole summation range (number of unit cell)
ex_range = 1; % Exchange interaction range (number of unit cell)
gLande_Ho = 1.25;
ELEf = gLande_Ho*muB*J2meV;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173
% NUCf = 4.732 * 3.1519e-5;   % Original values from linear response code
% NUCf = 1.668 * 3.15245e-5;  % http://easyspin.org/documentation/isotopetable.html
theta = 0; % angle between the external field and c-axis.

% freq_total = (0:0.5:160);
freq_total = (0:0.025:7);
for ii = 1:length(Temperatures)
    location = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I';
%     location = '/Volumes/GoogleDrive/My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I';
    filename = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg', Temperatures(ii), theta, phi),'.mat');
    file = fullfile(location,filename);
    load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" (eigenvalues) and "V" (eigenstates)
    fields = vecnorm(fff,2,1);
%     fields = fields(fields <= 5); % Truncate the field range
    fprintf('Calculating for T = %.3f K.\n', Temperatures(ii));
    %         [fields, freq_total, rechi, imchi] = linear_response(eee,fff,ttt,vvv);
    if Options.RPA == true
        switch Options.scanMode %1. Field plot with RPA. 2. wavevector plot with RPA
            case 'field'
                qz = 0.0;
%                 qz = 7.305; % Wavelength = 0.1369 m, frequency = 2.19 GHz
%                 qz = 12.175; % Wavelength = 0.0821 m, frequency = 3.65 GHz
%                 qx = [0.01 0.1 0.3 0.6 1]';
                qy = zeros(size(qz,1),1);
                qx = zeros(size(qz,1),1);
                qvec = [qx qy qz];
                eigenW = vvv;
                eigenE = eee;
            case 'qvec'
%                 qz = (6.9:0.0025:7.5)';
                qz = (0:0.0001:0.05)';
                qx = zeros(size(qz,1),1);
                qy = zeros(size(qz,1),1);
                qvec = [qx qy qz];
%                 B0 = 4.25; % critical field at 300 mK
                B0 = 5.25;
                bidx = int16.empty(0,length(B0));
                for jj = 1:length(B0)
                    [~,bidx(jj)] = min(abs(vecnorm(fields,1,1)-B0(jj)));
                end
                fields = fff(:,bidx);
                eigenW = vvv(bidx,:,:);
                eigenE = eee(bidx,:);
            otherwise
                disp('Unknown Scan Mode!')
                return
        end
        [fields, freq_total, chi0, ~] = linear_response(eigenE,fields,freq_total,ttt,eigenW,gama);
        [fields, freq_total, chiq.J] = RPA(qvec, fields, freq_total, ion, chi0.J); % Electronic susceptibilities
        [fields, freq_total, chiq.IJ] = RPA(qvec, fields, freq_total, ion, chi0.IJ); % cross term
        [fields, freq_total, chiq.I] = RPA(qvec, fields, freq_total, ion, chi0.I); % Nuclear susceptibilities

        chiq = ELEf^2*chiq.J + 2*ELEf*NUCf*chiq.IJ + NUCf^2*chiq.I; % Full electronuclear susceptibility expansion
%         chiq = ELEf^2*chiq.J + 2*ELEf*NUCf*chiq.IJ; % Electron susceptibility + the cross term
%         chiq = ELEf^2*chiq.J + NUCf^2*chiq.I; % Electron susceptibility + nuclear susceptibilities
%         chiq = ELEf^2*chiq.J; % Electron susceptibility only
%         chiq = NUCf^2*chiq.I; % Nuclear susceptibility only
    else
        eigenE = eee;
        eigenW = vvv;
        [fields, freq_total, chi0, ~] = linear_response(eigenE,fields,freq_total,ttt,eigenW,gama);
    end
    chi0 = ELEf^2*chi0.J + 2*ELEf*NUCf*chi0.IJ + NUCf^2*chi0.I; % Full electronuclear susceptibility expansion
%     chi0 = ELEf^2*chi0.J + 2*ELEf*NUCf*chi0.IJ; % Electron susceptibility + the cross term
%     chi0 = ELEf^2*chi0.J + NUCf^2*chi0.I; % Electron susceptibility + nuclear susceptibilities
%     chi0 = ELEf^2*chi0.J; % Electron susceptibility only
%     chi0 = NUCf^2*chi0.I; % Nuclear susceptibility only
    
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
        switch Options.scanMode
            case 'field'
                figs(fields,freq_total,chi0,gama,[0,0,0],Options,Temperatures(ii),"\chi_0");
                if Options.RPA == true
                    figs(fields,freq_total,chiq,gama,qvec,Options,Temperatures(ii),"\chi_{RPA}");
                    figs(fields,freq_total,chiq-chi0,gama,qvec,Options,Temperatures(ii),"\chi_{RPA}-\chi_0");
                end
            case 'qvec'
                chi0_p = repmat(chi0,1,1,1,length(qvec));
                chi0_p = permute(chi0_p,[1 2 3 5 4]); % Add the additional dimension for the plotting purpose
                figs(fields,freq_total,chi0_p,gama,qvec,Options,Temperatures(ii),"chi_0");
                figs(fields,freq_total,chiq,gama,qvec,Options,Temperatures(ii),"\chi_{RPA}");
        end
    end
end
end

function figs(fields,freq_total,chi,gama,qvec,Options,temperature,fig_tit)
global f2E
if Options.meV == true
    freq_total = freq_total*f2E;
end
pos0 = [100 300 600 400]; % initial figure position
pos_inc = [150 0 0 0];
indx = [strcat(fig_tit,"^{xx}") strcat(fig_tit,"^{yy}") strcat(fig_tit,"^{zz}")]; % index for figure legend
switch Options.scanMode
    case 'field'
        for ii = 1:size(qvec,1)
            % Color plot the imaginary part of the susceptibilities
            for jj = 3 % plot only chi_zz component
%             for jj = 1:3
                fig2 = figure;
                set(fig2,'position',pos0 + (jj-1)*pos_inc);
%                 hp2 = pcolor(fields(1,:),freq_total,squeeze(imag(chi(jj,jj,:,:,ii))));
                hp2 = pcolor(fields(1,:),freq_total,mag2db(abs(squeeze(imag(chi(jj,jj,:,:,ii))))));
                set(hp2, 'edgeColor','none')
                caxis([0 5]);
                colorbar
                legend(['\gamma =' num2str(gama,'%.2e meV')]);
                xlabel('Magnetic field (T)')
                ylabel('Frequency (GHz)')
                if Options.meV == true
                    ylabel('Energy (meV)')
                end
                title(['Im[', char(indx(jj)), '] at ', sprintf('T = %.3f K and ', temperature),sprintf('Q = [%1$.2f %2$.2f %3$.2f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])
                
                % Plot the real part of the susceptibility of the z component
                fig3 = figure;
                set(fig3,'position',pos0 + jj*pos_inc);
%                 hp3 = pcolor(fields(1,:),freq_total,squeeze(real(chi(jj,jj,:,:,ii))));
                hp3 = pcolor(fields(1,:),freq_total,mag2db(abs(squeeze(real(chi(jj,jj,:,:,ii))))));
                set(hp3, 'edgeColor','none')
                caxis([0 5]);
                colorbar
                legend(['\gamma =' num2str(gama,'%.2e meV')]);
                xlabel('Magnetic field (T)')
                ylabel('Frequency (GHz)')
                if Options.meV == true
                    ylabel('Energy (meV)')
                end
                title(['Re[', char(indx(jj)), '] at ', sprintf('T = %.3f K and ', temperature),sprintf('Q = [%1$.2f %2$.2f %3$.2f]',qvec(ii,1),qvec(ii,2),qvec(ii,3))])
            end
        end
    case 'qvec'
        qvec = qvec(:,3); % use qz
        for ii = 1:size(fields,2)
            for jj = 3 % plot only chi_zz component
%           for jj = 1:3
                % Color plot the imaginary part of the susceptibilities of z component
                fig2 = figure;
                set(fig2,'position',pos0 + (jj-1)*pos_inc);
                hp2 = pcolor(qvec,freq_total,real(squeeze(chi(jj,jj,:,ii,:))));
                set(hp2, 'edgeColor','none')
                set(gca, 'xdir', 'reverse' )
                caxis([0 100]);
                colorbar
                legend(['\gamma =' num2str(gama,'%.2e meV')]);
                xlabel(sprintf('Q = [h, 0, 0]'))
                ylabel('Frequency (GHz)')
                if Options.meV == true
                    ylabel('Energy (meV)')
                end
                title(['Re[', char(indx(jj)), '] at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])
                
                % Plot the real part of the susceptibility of the z component
                fig3 = figure;
                set(fig3,'position',pos0 + jj*pos_inc);
                hp3 = pcolor(qvec,freq_total,imag(squeeze(chi(jj,jj,:,ii,:))));
                set(hp3, 'edgeColor','none')
                set(gca, 'xdir', 'reverse' )
                caxis([0 100]);
                colorbar
                legend(['\gamma =' num2str(gama,'%.2e meV')]);
                xlabel(sprintf('Q = [h, 0, 0]'))
                ylabel('Frequency (GHz)')
                if Options.meV == true
                    ylabel('Energy (meV)')
                end
                title(['Im[', char(indx(jj)), '] at ',sprintf('T = %.3f K. and ', temperature),sprintf('B = [%1$.2f %2$.2f %3$.2f] T.',fields(1,ii),fields(2,ii),fields(3,ii))])
            end
        end
    otherwise
        disp("Unknow scan mode!")
        return
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

function [fields, freq_total, chi0, JIz_exp] = linear_response(eigenE,fields,freq_total,temperature,eigenW,gama)
global ELEf NUCf f2E
% Declare susceptibility tensor
chi0r_J = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0i_J = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0r_I = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0i_I = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0r_IJ = zeros(3,3,length(freq_total(1,:)),size(fields,2));
chi0i_IJ = zeros(3,3,length(freq_total(1,:)),size(fields,2));

%Initiate J operators
J = 8; % Ho
I = 3.5; % Ho
Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
JhT.z=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*I+1)); % Expand to match the dimension of Hilbert space
Jmh=kron(Jm,eye(2*I+1));
JhT.x=(Jph+Jmh)/2;
JhT.y=(Jph-Jmh)/2i;

%Initiate I operators
Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
IhT.z=kron(eye(2*J+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*J+1),Ip); % Expand to match the dimension of Hilbert space
Imh=kron(eye(2*J+1),Im);
IhT.x=(Iph+Imh)/2;
IhT.y=(Iph-Imh)/2i;
        
% IJ_hT.x = JhT.x+NUCf/ELEf*IhT.x; % Hybridized electronuclear spin operator
% IJ_hT.y = JhT.y+NUCf/ELEf*IhT.y;
% IJ_hT.z = JhT.z+NUCf/ELEf*IhT.z;

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

for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (m);
    omega = freq*f2E;   % define frequency sweep range (meV)
%     for k = 1:size(fields,2) % for debugging: calculate susceptibility for all fields
    parfor k = 1:size(fields,2) % calculate susceptibility for all fields
        v = squeeze(eigenW(k,:,:)); % Obtain the corresponding eigen vectors
        en = squeeze(eigenE(k,:)); % Obtain the corresponding eigen energies [meV]
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
        deno1 = gamma ./ (EE.^2 + gamma.^2);
        deno2 = EE ./ (EE.^2 + gamma.^2);

%         [chi0r_J(:,:,m,k), chi0i_J(:,:,m,k)] = chi_Mx(IJ_hT,IJ_hT,v,NN,deno1,deno2); % hyperdized electronuclear operator

        [chi0r_J(:,:,m,k), chi0i_J(:,:,m,k)] = chi_Mx(JhT,JhT,v,NN,deno1,deno2); % Electornic spin operators
        [chi0r_IJ(:,:,m,k), chi0i_IJ(:,:,m,k)] = chi_Mx(IhT,JhT,v,NN,deno1,deno2); % Electro-Nuclear cross term
        [chi0r_I(:,:,m,k), chi0i_I(:,:,m,k)] = chi_Mx(IhT,IhT,v,NN,deno1,deno2); % Nuclear spin operators
    end
end

chi0.J = (chi0r_J + 1i*chi0i_J); % Electronic susceptibilities
chi0.I = (chi0r_I + 1i*chi0i_I); % Nuclear susceptibilities
chi0.IJ = (chi0r_IJ + 1i*chi0i_IJ); % Electronuclear Cross term
end

function [xr,xi] = chi_Mx(opr1, opr2, wav, NN, deno1, deno2)
% Calculation of matrix elements for the real and imaginary susceptibilities
% opr1: operator 1
% opr2: operator 2
% wav: eigenfunction of the hamiltonian
% NN: population differene between two eigenstates
% deno1: denominator for real part of the susceptibility
% deno2: denominator for imaginary part of the susceptibility

ttx1  = wav'  * opr1.x * wav;
tty1  = wav'  * opr1.y * wav;
ttz1  = wav'  * opr1.z * wav;

ttx2  = wav'  * opr2.x * wav;
tty2  = wav'  * opr2.y * wav;
ttz2  = wav'  * opr2.z * wav;

% Calculate susceptibilities along a-axis
chi_tx  = (ttx1) .* (ttx2.') .* NN .* deno1; % Imaginary part
chi_t1x = (ttx1) .* (ttx2.') .* NN .* deno2; % Real part
xx = sum(sum(chi_tx));
xx1 = sum(sum(chi_t1x));

% Calculate susceptibilities along b-axis
chi_ty  = (tty1) .* (tty2.') .* NN .* deno1;
chi_t1y = (tty1) .* (tty2.') .* NN .* deno2;
yy = sum(sum(chi_ty));
yy1 = sum(sum(chi_t1y));

% Calculate susceptibilities along c-axis
chi_tz  = (ttz1) .* (ttz2.') .* NN .* deno1;
chi_t1z = (ttz1) .* (ttz2.') .* NN .* deno2;
zz = sum(sum(chi_tz));
zz1 = sum(sum(chi_t1z));

% Calculate susceptibilities along ab-axis
chi_txy  = (ttx1) .* (tty2.') .* NN .* deno1;
chi_t1xy = (ttx1) .* (tty2.') .* NN .* deno2;
xy = sum(sum(chi_txy));
xy1 = sum(sum(chi_t1xy));

% Calculate susceptibilities along ac-axis
chi_txz  = (ttx1) .* (ttz2.') .* NN .* deno1;
chi_t1xz = (ttx1) .* (ttz2.') .* NN .* deno2;
xz = sum(sum(chi_txz));
xz1 = sum(sum(chi_t1xz));

% Calculate susceptibilities along ba-axis
chi_tyx  = (tty1) .* (ttx2.') .* NN .* deno1;
chi_t1yx = (tty1) .* (ttx2.') .* NN .* deno2;
yx = sum(sum(chi_tyx));
yx1 = sum(sum(chi_t1yx));

% Calculate susceptibilities along bc-axis
chi_tyz  = (tty1) .* (ttz2.') .* NN .* deno1;
chi_t1yz = (tty1) .* (ttz2.') .* NN .* deno2;
yz = sum(sum(chi_tyz));
yz1 = sum(sum(chi_t1yz));

% Calculate susceptibilities along ca-axis
chi_tzx  = (ttz1) .* (ttx2.') .* NN .* deno1;
chi_t1zx = (ttz1) .* (ttx2.') .* NN .* deno2;
zx = sum(sum(chi_tzx));
zx1 = sum(sum(chi_t1zx));

% Calculate susceptibilities along cb-axis
chi_tzy  = (ttz1) .* (tty2.') .* NN .* deno1;
chi_t1zy = (ttz1) .* (tty2.') .* NN .* deno2;
zy = sum(sum(chi_tzy));
zy1 = sum(sum(chi_t1zy));

xr = [xx1 xy1 xz1
      yx1 yy1 yz1
      zx1 zy1 zz1]; % Real part of susceptibility

xi = [xx xy xz
      yx yy yz
      zx zy zz]; % Imaginary part of susceptibility
end

function [fields, freq_total, chiq] = RPA(qvec, fields, freq_total, ion, chi0)

global muB mu0 J2meV gLande_Ho dip_range ex_range
N = 4; % Number of magnetic atoms in unit cell
a = [5.175 0 0;
     0 5.175 0;
     0 0 10.75]; % Lattice constant for LiHoF4
gfac = (muB)^2*(mu0/4/pi)*J2meV*10^30;

chiq = zeros(3,3,length(freq_total(1,:)),size(fields,2),size(qvec,1));
D = zeros(3,3,N,N,size(qvec,1));
for jj=1:size(qvec,1)
    D(:,:,:,:,jj) = gLande_Ho^2*(gfac*dipole_direct(qvec(jj,:),dip_range,a))...
        + exchange(qvec(jj,:),ion.ex(2),a,ex_range);
end

deno = zeros(3,3,size(freq_total,1),size(fields,2)); % RPA correction factor (denominator)
for ii = 1:size(fields,2)
    for nq = 1:size(qvec,1)
        Jq = sum(sum(D(:,:,:,:,nq),4),3)/4; % average over the four sites in the unit cell and avoid double counting
        parfor kk = 1:length(freq_total(1,:))
            MM = chi0(:,:,kk,ii).*Jq;
            deno(:,:,kk,ii) = (ones(size(MM))- MM); %#ok<PFOUS> Suppress parfor warning for this line
            chiq(:,:,kk,ii,nq) = squeeze(deno(:,:,kk,ii)).\chi0(:,:,kk,ii);
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