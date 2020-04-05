function linear_response_theory_publish
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output')
clearvars;
% Temperatures = [0.080, 0.160, 0.200, 0.300, 0.500];
Temperatures=[0.150];
filenames = strcat(num2str(Temperatures(:),'%.3f'),'.mat');
load(filenames,'-mat','eee','fff','ttt','vvv'); % loads variables "fields", "temp", "E" and "V" 
% which are eigenstates and eigenvalues calculated in the mean-field model 
% as a function of transverse field and temperature

E = eee;
V = vvv;
temp = ttt;
fields = vecnorm(fff);
freq_total = (1:0.02:5);

imchix = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);
imchiy = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);
imchiz = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);

rechi1x = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);
rechi1y = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);
rechi1z = double.empty(length(temp(1,:)),length(freq_total(1,:)),length(fields(1,:)),0);

for l = 1:length(temp(1,:)) % calculate susceptibility for all temperatures
    tt = temp(1,l);

    for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
        freq = freq_total (m);
        J=8;
        I=3.5;
        gLande_Ho=1.25;

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
        Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1); % Nuclear spin ladder operator
        Im=Ip'; % Nuclear spin ladder operator
        Iph=kron(eye(2*J+1),Ip); % Expand to match the Hilbert space
        Imh=kron(eye(2*J+1),Im); % Expand to match the Hilbert space
        Ixh=(Iph+Imh)/2;
        Iyh=(Iph-Imh)/2i;

        ghztomeV = 1/241.8;
        omega = freq*ghztomeV;     % define frequency sweep range (meV)
        gama = 0.00005; % define lifetime (meV) (theoretical basis?)

        for k = 1:length(fields(1,:)) % calculate susceptibility for all fields
            v = squeeze ( V(k,l,:,:) ); % Obtain the corresponding eigen vectors
            en = squeeze ( E(k,l,:) ); % Obtain the corresponding eigen energies in meV
%             field = fields(1,k); % Obtain the corresponding field

%             N = length(en);
%             chi_t = zeros(1,N^2);
            %ll = 1;
%             zz = double.empty(0,N);
%             beta = 1/(tt/11.6);
            beta = 1/(tt*8.617E-2); %[meV^-1]
            z=sum(exp(-beta*en));
            zz=exp(-beta*en)/z;
            [n,np]=meshgrid(zz,zz);
            NN=n-np;
            [ee,eep]=meshgrid(en,en);
%             EE1=1./(ee-eep-omega);
            EE = eep-ee-omega;
            gamma = ones(size(EE))*gama;
            G = gamma ./ (EE.^2 + gamma.^2);
            G1 = EE ./ (EE.^2 + gamma.^2);

            ELEf = gLande_Ho * 0.05788;     % Lande factor x Bohr magneton (meV T^-1)
            NUCf = 4.732 * 3.1519e-5;   % Nuclear Lande factor
            JxhT = Jxh * ELEf;
            IxhT = Ixh * NUCf;
            
            JyhT = Jyh * ELEf;
            IyhT = Iyh * NUCf;
            
            JzhT = Jzh * ELEf;
            IzhT = Izh * NUCf;

% Calculate susceptibilities along x-axis
            ttx  = v'  * (JxhT+IxhT) * v; 
            chi_tx  = (ttx) .* (ttx.') .* NN .* G;
            chi_t1x = (ttx) .* (ttx.') .* NN .* G1;
            sss=sum(sum(chi_tx)); 
            sss1=sum(sum(chi_t1x));
            imchix  (l,m,k,1) =  real(sss)   ;
            rechi1x (l,m,k,1) =  real(sss1)  ; 

% Calculate susceptibilities along y-axis
            tty  = v'  * (JyhT+IyhT) * v; 
            chi_ty  = (tty) .* (tty.') .* NN .* G;
            chi_t1y = (tty) .* (tty.') .* NN .* G1;
            sss=sum(sum(chi_ty));
            sss1=sum(sum(chi_t1y));
            imchiy  (l,m,k,1) =  real(sss)   ;
            rechi1y (l,m,k,1) =  real(sss1)  ;   
%             titttz = 'S(Jyy+Iyy)';

% Calculate susceptibilities along c-axis
            ttz  = v'  * (JzhT+IzhT) * v;
            chi_tz  = (ttz) .* (ttz.') .* NN .* G;
            chi_t1z = (ttz) .* (ttz.') .* NN .* G1;
            sss=sum(sum(chi_tz));
            sss1=sum(sum(chi_t1z));
            imchiz  (l,m,k,1) =  real(sss)   ;
            rechi1z (l,m,k,1) =  real(sss1)  ;
            %titttz = 'S(Jzz+Izz)';
        end
    end
end

%% Save the susceptibilities
x1x = squeeze(rechi1x); 
x2x = squeeze(imchix);
save(strcat('x1x_x2x_',num2str(Temperatures(:)*1000),'mK'),'fields','freq_total','x1x','x2x')

x1y = squeeze(rechi1y); 
x2y = squeeze(imchiy);
save(strcat('x1y_x2y_',num2str(Temperatures(:)*1000),'mK'),'fields','freq_total','x1y','x2y')

x1z = squeeze(rechi1z); 
x2z = squeeze(imchiz);
save(strcat('x1z_x2z_',num2str(Temperatures(:)*1000),'mK'),'fields','freq_total','x1z','x2z')
%% Color plot the susceptibilities
% % Color plot of the imaginary part of the susceptibility of x component
%     hfig0 = figure (1);
%     clf
%     hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchix)));
%     set(hp1, 'edgeColor','none')
%     caxis([-23 2]);
%     colorbar
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     title({'Imaginary part of Susceptibility (log scale)', 'in x direction'})
%     
% % Color plot of the imaginary part of the susceptibility of y component
%     hfig1 = figure (2);
%     clf
%     hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchiy)));
%     set(hp1, 'edgeColor','none')
%     caxis([-23 2]);
%     colorbar
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     title({'Imaginary part of Susceptibility (log scale)', 'in y direction'})
%     
% % Color plot the imaginary part of the susceptibilities of z component
%     hfig2 = figure (3);
%     clf
%     hp2 = pcolor(fields(1,:),freq_total,squeeze(log(imchiz)));
%     set(hp2, 'edgeColor','none')
%     caxis([-23 2]);
%     colorbar
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     title({'Imaginary part of Susceptibility (log scale)', 'in z direction'})
% 
% % Plot the real part of the susceptibility of the z component
%     hfig4 = figure (4);
%     clf
%     set(hfig1,'position',[50 100 600 400])
%     hp2 = pcolor(fields(1,:),freq_total,squeeze(rechi1z));
%     set(hp2, 'edgeColor','none')
%     caxis([-23 2]);
%     colorbar
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     title({'Real part of Susceptibility', 'in z direction'})
end