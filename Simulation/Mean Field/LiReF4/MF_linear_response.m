function MF_linear_response
clearvars;
Options.RPA = false; % Apply random phase approximation (RPA) correction
Options.plotting = false; % Decide whether or not to plot the data at the end
Options.saving = true;

% Temperatures = [0.08 0.1 0.12 0.15 0.2 0.24];
Temperatures = 0.08;
theta = 0.0;
phi = 0.0;
% gama = 0.0001; % define lifetime (meV)
gama = 1e-3; % hyperfine state lifetime (meV)

    for ii = 1:length(Temperatures)
        location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output\without Hz_I';
        filename = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg', Temperatures(ii), theta, phi),'.mat');
        file = fullfile(location,filename);
        load(file,'-mat','eee','fff','ttt','vvv','ion'); % loads variables "fields", "temp", "E" and "V" 
        % which are eigenstates and eigenvalues calculated in the mean-field model 
        % as a function of transverse field and temperature
        
        fprintf('Calculating for T = %.3f.\n', Temperatures(ii));
%         [fields, freq_total, rechi, imchi] = linear_response(eee,fff,ttt,vvv);
        if Options.RPA == true
            dip_range = 80;
            qvec = [0 0 1];
            [fields, freq_total, rechi0, imchi0, chi0r, chi0i, ~] = linear_response(eee,fff,ttt,vvv,gama);
            [fields, freq_total, rechi, imchi] = RPA(qvec, fields, freq_total, ion, chi0r, chi0i, dip_range);
        else
            [fields, freq_total, rechi0, imchi0, ~, ~, ~] = linear_response(eee,fff,ttt,vvv,gama);
        end
        
        if Options.saving == true % Save the susceptibilities
            file_name1 = strcat('LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e.mat', ttt, theta, phi,gama));
            savefile1 = fullfile(location,file_name1);
            if Options.RPA == true                
                file_name2 = strcat('RPA_LiHoF4_x1z_x2z_',sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg_%4$.2e.mat', ttt, theta, phi,gama));
                savefile2 = fullfile(location,file_name2);
                save_file(fields,freq_total,rechi0,imchi0,gama,savefile1); % Save the data without RPA corrections
                save_file(fields,freq_total,rechi1,imchi,gama,savefile2); % save the data with RPA corrections
            else
                save_file(fields,freq_total,rechi0,imchi0,gama,savefile1); % Save the data without RPA corrections
            end
        end
        
        if Options.plotting == true % Plot the susceptibilities
            if Options.RPA == true
                plot(fields,freq_total,rechi0,imchi0,gama); % Plot the data without RPA corrections
                plot(fields,freq_total,rechi,imchi,gama); % Plot the data with RPA corrections
            else
                plot(fields,freq_total,rechi0,imchi0,gama); % Plot the data without RPA corrections
            end
        end
    end
end

function plot(fields,freq_total,rechi,imchi,gama)
pos0 = [100 300 600 400]; % initial figure position
pos_inc = [100 0 0 0];
% Color plot of the imaginary part of the susceptibility of x component
fig0 = figure;
set(fig0,'position',pos0);
hp0 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.x)));
set(hp0, 'edgeColor','none')
caxis([-23 2]);
colorbar
legend(['$\gamma$ =' num2str(gama,'%.2e')]);
xlabel('Magnetic field (T)')
ylabel('Frequency (GHz)')
title({'Imaginary part of $\chi$ (log scale) in x direction'})

% Color plot of the imaginary part of the susceptibility of y component
fig1 = figure;
set(fig1,'position',pos0 + pos_inc);
hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.y)));
set(hp1, 'edgeColor','none')
caxis([-23 2]);
colorbar
legend(['$\gamma$ =' num2str(gama,'%.2e')]);
xlabel('Magnetic field (T)')
ylabel('Frequency (GHz)')
title({'Imaginary part of $\chi$ (log scale) in y direction'})

% Color plot the imaginary part of the susceptibilities of z component
fig2 = figure;
set(fig2,'position',pos0 + 2*pos_inc);
hp2 = pcolor(fields(1,:),freq_total,squeeze(log(imchi.z)));
set(hp2, 'edgeColor','none')
caxis([-23 2]);
colorbar
legend(['$\gamma$ =' num2str(gama,'%.2e')]);
xlabel('Magnetic field (T)')
ylabel('Frequency (GHz)')
title({'Imaginary part of $\chi$ (log scale) in z direction'})

% Plot the real part of the susceptibility of the z component
fig3 = figure;
set(fig3,'position',pos0 + 3*pos_inc);
hp3 = pcolor(fields(1,:),freq_total,squeeze(rechi.z));
set(hp3, 'edgeColor','none')
caxis([-23 2]);
colorbar
legend(['$\gamma$ =' num2str(gama,'%.2e')]);
xlabel('Magnetic field (T)')
ylabel('Frequency (GHz)')
title({'Real part of $\chi$ in z direction'})

% % Plot the expectation value of <Jz+Iz>
%     fig4 = figure;
%     set(fig4,'position',pos0 + 4*pos_inc);
%     hp4 = plot(fields(1,:),squeeze(Jz_exp),'-o');
%     xlabel('Magnetic field (T)')
%     ylabel('<J_z + I_z>')
%     title({'Expectation value of J_z + I_z', 'in z direction'})
end
function save_file(fields,freq_total,rechi,imchi,gama,savefile)
% x1x = squeeze(rechi1.x);
% x2x = squeeze(imchi.x);
% save(strcat('LiHoF4_x1x_x2x_',sprintf('%1$3.3fK_%2$uDeg', ttt, theta),'fields','freq_total','x1x','x2x'));
%
% x1y = squeeze(rechi1.y);
% x2y = squeeze(imchi.y);
% save(strcat('LiHoF4_x1y_x2y_',sprintf('%1$3.3fK_%2$uDeg', ttt, theta),'fields','freq_total','x1y','x2y'));

x1z = squeeze(rechi.z);
x2z = squeeze(imchi.z);
save(savefile,'fields','freq_total','x1z','x2z','gama','-v7.3');
end
function [fields, freq_total, rechi0, imchi0, chi0r, chi0i, JIz_exp]=linear_response(eee,fff,ttt,vvv,gama)
% Calculation of susceptibilities
E = eee;
V = vvv;
fields = vecnorm(fff);
freq_total = (1:0.01:5);
% freq_total = (1:0.01:5);

imchix = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiy = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiz = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

rechi1x = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1y = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1z = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

chi0r = zeros(3,3,length(freq_total(1,:)),length(fields(1,:)));
chi0i = zeros(3,3,length(freq_total(1,:)),length(fields(1,:)));

J=8;
I=3.5;
gLande_Ho=1.25;
ELEf = gLande_Ho * 0.05788;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173
% NUCf = 4.732 * 3.1519e-5;   % Original code
    
%Initiate J operators
Jz=diag(J:-1:-J); % Jz = -J, -J+1,...,J-1,J
Jzh=kron(Jz,eye(2*I+1)); % Expand Jz space to include nuclear degree of freedom
Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1); % electronic spin ladder operator
Jm=Jp'; % electronic spin ladder operator
Jph=kron(Jp,eye(2*I+1)); % Expand Hilbert space
Jmh=kron(Jm,eye(2*I+1)); % Expand Hilbert space
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;
%tensor product of cristal field to include nuclear moments
%Initiate I operators
Iz=diag(I:-1:-I); %Iz = -I, -I+1,...,I-1,I
Izh=kron(eye(2*J+1),Iz); % Expand Hilbert space
Ip=diag(sqrt((I-((I-1):-1:-I)).*(I+1+((I-1):-1:-I))),1); % Nuclear spin ladder operator
Im=Ip'; % Nuclear spin ladder operator
Iph=kron(eye(2*J+1),Ip); % Expand to match the Hilbert space
Imh=kron(eye(2*J+1),Im); % Expand to match the Hilbert space
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

% Single out <Jz+Iz> calculations
% Jz_exp = double.empty(0,length(fields(1,:)));
JIz_exp = double.empty(0,length(fields(1,:)));
for kk = 1:length(fields(1,:)) % calculate susceptibility for all fields
    v = squeeze(squeeze(V(kk,:,:,:))); % Obtain the corresponding eigen vectors
    JzhT = Jzh * ELEf;
    IzhT = Izh * NUCf;
%     tz = v'  * JzhT * v;
    ttz  = v'  * (JzhT+IzhT) * v;
%     Jz_exp(1,kk) = sqrt(sum(sum((tz) .* (tz.'))));
    JIz_exp(1,kk) = sqrt(sum(sum((ttz) .* (ttz.'))));
end

for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (m);
    f2E = 1/241.8;  % GHz to meV
    omega = freq*f2E;   % define frequency sweep range (meV)
%     for k = 1:length(fields(1,:)) % for debugging: calculate susceptibility for all fields
    parfor k = 1:length(fields(1,:)) % calculate susceptibility for all fields
        v = squeeze ( squeeze(V(k,:,:,:)) ); % Obtain the corresponding eigen vectors
        en = squeeze ( squeeze(E(k,:,:)) ); % Obtain the corresponding eigen energies in meV
%       N = length(en);
%       chi_t = zeros(1,N^2);
%       ll = 1;
        beta = 11.6/ttt; %[meV^-1]
        z=sum(exp(-beta*en));
        zz=exp(-beta*en)/z;
        [n,np]=meshgrid(zz,zz);
        NN=n-np;
        [ee,eep]=meshgrid(en,en);
        EE = eep-ee-omega;
        gamma = ones(size(EE))*gama;
        G = gamma ./ (EE.^2 + gamma.^2); 
        G1 = EE ./ (EE.^2 + gamma.^2);  

        JxhT = Jxh * ELEf;
        IxhT = Ixh * NUCf;

        JyhT = Jyh * ELEf;
        IyhT = Iyh * NUCf;

        JzhT = Jzh * ELEf;
        IzhT = Izh * NUCf;

% Calculate susceptibilities along a-axis
        ttx  = v'  * (JxhT+IxhT) * v; 
        chi_tx  = (ttx) .* (ttx.') .* NN .* G;
        chi_t1x = (ttx) .* (ttx.') .* NN .* G1;
        xx=sum(sum(chi_tx)); 
        xx1=sum(sum(chi_t1x));
        imchix  (m,k,1) =  real(xx)   ;
        rechi1x (m,k,1) =  real(xx1)  ; 

% Calculate susceptibilities along b-axis
        tty  = v'  * (JyhT+IyhT) * v; 
        chi_ty  = (tty) .* (tty.') .* NN .* G;
        chi_t1y = (tty) .* (tty.') .* NN .* G1;
        yy=sum(sum(chi_ty));
        yy1=sum(sum(chi_t1y));
        imchiy  (m,k,1) =  real(yy)   ;
        rechi1y (m,k,1) =  real(yy1)  ;   

% Calculate susceptibilities along c-axis
        ttz  = v'  * (JzhT+IzhT) * v;
        chi_tz  = (ttz) .* (ttz.') .* NN .* G;
        chi_t1z = (ttz) .* (ttz.') .* NN .* G1;
        zz=sum(sum(chi_tz));
        zz1=sum(sum(chi_t1z));
        imchiz  (m,k,1) =  real(zz)   ;
        rechi1z (m,k,1) =  real(zz1)  ;
        
% Calculate susceptibilities along ab-axis
        chi_txy  = (ttx) .* (tty.') .* NN .* G;
        chi_t1xy = (ttx) .* (tty.') .* NN .* G1;
        xy=sum(sum(chi_txy)); 
        xy1=sum(sum(chi_t1xy));
%         imchixy  (m,k,1) =  real(xy)   ;
%         rechi1xy (m,k,1) =  real(xy1)  ; 

% Calculate susceptibilities along ac-axis
        chi_txz  = (ttx) .* (ttz.') .* NN .* G;
        chi_t1xz = (ttx) .* (ttz.') .* NN .* G1;
        xz=sum(sum(chi_txz));
        xz1=sum(sum(chi_t1xz));
%         imchixz  (m,k,1) =  real(xz)   ;
%         rechi1xz (m,k,1) =  real(xz1)  ;   

% Calculate susceptibilities along ba-axis
        chi_tyx  = (tty) .* (ttx.') .* NN .* G;
        chi_t1yx = (tty) .* (ttx.') .* NN .* G1;
        yx=sum(sum(chi_tyx));
        yx1=sum(sum(chi_t1yx));
%         imchiyx  (m,k,1) =  real(yx)   ;
%         rechi1yx (m,k,1) =  real(yx1)  ;
        
% Calculate susceptibilities along bc-axis
        chi_tyz  = (tty) .* (ttz.') .* NN .* G;
        chi_t1yz = (tty) .* (ttz.') .* NN .* G1;
        yz=sum(sum(chi_tyz));
        yz1=sum(sum(chi_t1yz));
%         imchiyz  (m,k,1) =  real(yz)   ;
%         rechi1yz (m,k,1) =  real(yz1)  ;
        
% Calculate susceptibilities along ca-axis
        chi_tzx  = (ttz) .* (ttx.') .* NN .* G;
        chi_t1zx = (ttz) .* (ttx.') .* NN .* G1;
        zx=sum(sum(chi_tzx));
        zx1=sum(sum(chi_t1zx));
%         imchizx  (m,k,1) =  real(zx)   ;
%         rechi1zx (m,k,1) =  real(zx1)  ;
        
% Calculate susceptibilities along cb-axis
        chi_tzy  = (ttz) .* (tty.') .* NN .* G;
        chi_t1zy = (ttz) .* (tty.') .* NN .* G1;
        zy=sum(sum(chi_tzy));
        zy1=sum(sum(chi_t1zy));
%         imchizy  (m,k,1) =  real(zy)   ;
%         rechi1zy (m,k,1) =  real(zy1)  ;

        chi0i(:,:,m,k)=[xx xy xz
            yx yy yz
            zx zy zz];
        
        chi0r(:,:,m,k)=[xx1 xy1 xz1
            yx1 yy1 yz1
            zx1 zy1 zz1];
    end
end
imchi0.x = imchix;
imchi0.y = imchiy;
imchi0.z = imchiz;

rechi0.x = rechi1x;
rechi0.y = rechi1y;
rechi0.z = rechi1z;
end
function [fields, freq_total, rechi, imchi] = RPA(qvec, fields, freq_total, ion, chi0r, chi0i, dip_range)

N = 4; % Number of magnetic atoms in unit cell
a = [5.175 0 0; 
     0 5.175 0;
     0 0 10.75]; % Lattice constant for LiHoF4
 
chir = zeros(3,3,length(freq_total(1,:)),size(qvec,1),length(fields(1,:)));
chii = zeros(3,3,length(freq_total(1,:)),size(qvec,1),length(fields(1,:)));
D = zeros(3,3,N,N,size(qvec,1));

imchix = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiy = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
imchiz = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

rechi1x = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1y = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);
rechi1z = double.empty(length(freq_total(1,:)),length(fields(1,:)),0);

for jj=1:size(qvec,1)
     D(:,:,:,:,jj)=dipole_direct(qvec(jj,:),dip_range,a) + exchange(qvec(jj,:),ion.ex(2),a);
end
D=D*(6/5)^2*0.05368;

for ii = 1:length(fields(1,:))
    for kk = 1:length(freq_total(1,:))
        for nq = 1:size(qvec,1)
            MM = zeros(3*N);
            NN = zeros(3*N);
            for n = 1:N
                for m = 1:N
                    MM((n-1)*3+(1:3),(m-1)*3+(1:3)) = chi0r(:,:,kk,ii)*D(:,:,n,m,nq);
                    NN((n-1)*3+(1:3),(m-1)*3+(1:3)) = chi0i(:,:,kk,ii)*D(:,:,n,m,nq);
                end
            end
            
            chir(:,:,kk,nq,ii) = 1/4*[eye(3) eye(3) eye(3) eye(3)]*...
                ((eye(size(MM))-MM)\([chi0r(:,:,kk,ii);chi0r(:,:,kk,ii);chi0r(:,:,kk,ii);chi0r(:,:,kk,ii)]));
            chii(:,:,kk,nq,ii) = 1/4*[eye(3) eye(3) eye(3) eye(3)]*...
                ((eye(size(NN))-NN)\([chi0i(:,:,kk,ii);chi0i(:,:,kk,ii);chi0i(:,:,kk,ii);chi0i(:,:,kk,ii)]));
            
            rechi1x(kk,ii,1) =  real(chir(1,1,kk,nq,ii));
            imchix(kk,ii,1) =  real(chii(1,1,kk,nq,ii));
            
            rechi1y(kk,ii,1) =  real(chir(2,2,kk,nq,ii));
            imchiy(kk,ii,1) =  real(chii(2,2,kk,nq,ii));
            
            rechi1z(kk,ii,1) =  real(chir(3,3,kk,nq,ii));
            imchiz(kk,ii,1) =  real(chii(3,3,kk,nq,ii));
        end
    end
end
rechi.x = rechi1x;
rechi.y = rechi1y;
rechi.z = rechi1z;

imchi.x = imchix;
imchi.y = imchiy;
imchi.z = imchiz;
end
