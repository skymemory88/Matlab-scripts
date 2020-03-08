function linear_response_theory
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output')
% Temperatures = [0.080, 0.160, 0.200, 0.300, 0.500];
% filenames = num2str(Temperatures,'.mat')
load('0.085.mat') % loads variables "fields", "temp", "E" and "V" 
% which are eigenstates and eigenvalues calculated in the mean-field model 
% as a function of transverse field and temperature

E = eee;
V = vvv;
temp = ttt;
fields = vecnorm(fff);
freq_total = (1:0.05:5);

for l = 1:length(temp(1,:)) % calculate susceptibility for all temperatures
    t = temp(1,l);

    for m = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
        freq = freq_total (m);
        J=8;
        I=3.5;
        % gLande_Ho=1.25;

        %Initiate J operators
        Jz=diag(J:-1:-J);
        Jzh=kron(Jz,eye(2*I+1));
        Jp=diag(sqrt((J-((J-1):-1:-J) ).*(J+1+( (J-1):-1:-J) )),1);
        Jm=Jp';
        Jph=kron(Jp,eye(2*I+1));
        Jmh=kron(Jm,eye(2*I+1));
        Jxh=(Jph+Jmh)/2;
        Jyh=(Jph-Jmh)/2i;
        %tensor product of cristal field to include nuclear moments
        %Initiate I operators
        Iz=diag(I:-1:-I);
        Izh=kron(eye(2*J+1),Iz);
        Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
        Im=Ip';
        Iph=kron(eye(2*J+1),Ip);
        Imh=kron(eye(2*J+1),Im);
        Ixh=(Iph+Imh)/2;
        Iyh=(Iph-Imh)/2i;

        ghztomeV = 1/241.8;
        omega = freq*ghztomeV;     % define frequency sweep range (meV)
        gama = 0.0001; % define lifetime (meV)

        for k = 1:length(fields(1,:)) % calculate susceptibility for all fields
            v = squeeze ( V(k,l,:,:) ); % Obtain the corresponding eigen vectors
            en = squeeze ( E(k,l,:) ); % Obtain the corresponding eigen energies
            field = fields(1,k); % Obtain the corresponding field

            N = length(en);
            chi_t = zeros(1,N^2);
%             ll = 1;
            zz = zeros(1,N);
            beta = 1/(t/11.6);
            z=sum(exp(-beta*en));
            zz=exp(-beta*en)/z;
            [n,np]=meshgrid(zz,zz);
            NN=n-np;
            [ee,eep]=meshgrid(en,en);
            EE1=1./(ee-eep-omega);
            EE = eep-ee-omega;
            gamma = ones(size(EE))*gama;
            G = gamma ./ (EE.^2 + gamma.^2);
            G1 = EE ./ (EE.^2 + gamma.^2);

            ELEf = 1.250 * 0.05788;
            NUCf = 4.732 * 3.1519e-5;
            JxhT = Jxh * ELEf;
            IxhT = Ixh * NUCf;
            JyhT = Jyh * ELEf;
            IyhT = Iyh * NUCf;
            JzhT = Jzh * ELEf;
            IzhT = Izh * NUCf;
            
%           Calculate susceptibilities along a(b)-axis
            tty  = v'  * (JyhT+IyhT) * v; 
            chi_ty  = (tty) .* (tty.') .* NN .* G;
            chi_t1y = (tty) .* (tty.') .* NN .* G1;
            sss=sum(sum(chi_ty));
            sss1=sum(sum(chi_t1y));
            imchiy  (l,m,k) =  real(sss)   ;
            rechi1y (l,m,k) =  real(sss1)  ;   
%             titttz = 'S(Jyy+Iyy)';

%           Calculate susceptibilities along c-axis
            ttz  = v'  * (JyhT+IyhT) * v;
            chi_tz  = (ttz) .* (ttz.') .* NN .* G;
            chi_t1z = (ttz) .* (ttz.') .* NN .* G1;
            sss=sum(sum(chi_tz));
            sss1=sum(sum(chi_t1z));
            imchiz  (l,m,k) =  real(sss)   ;
            rechi1z (l,m,k) =  real(sss1)  ;
%             titttz = 'S(Jzz+Izz)';

        end
        %         hfig1 = figure (1);
        %         clf
        %         set(hfig1,'position',[50 100 600 400])
        %         h1=plot (fields(1,:), imchi ,'r','LineWidth',2); 
        %         set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        %         set(gca,'fontsize',15)
        %         xlim([0 9]);
        %         xlabel('Magnetic field (T)','FontSize',15)
        %         ylabel('\chi'''' (arb. u.)','FontSize',15)
        % 
        %         hfig2 = figure (2);
        %         clf
        %         set(hfig2,'position',[680 100 600 400])
        %         h2=plot (fields(1,:), rechi1 ,'r','LineWidth',2); 
        %         set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        %         set(gca,'fontsize',15)
        %         xlim([0 9]);
        %         xlabel('Magnetic field (T)','FontSize',15)
        %         ylabel('\chi'' (arb. u.)','FontSize',15)
    end
end

% % Save the susceptibilities
x1y = squeeze(rechi1y); 
x2y = squeeze(imchiy);
save('x1y_x2y_85mK','fields','freq_total','x1y','x2y')

x1z = squeeze(rechi1z); 
x2z = squeeze(imchiz);
save('x1z_x2z_85mK','fields','freq_total','x1z','x2z')

% Color plot the susceptibilities
    hfig1 = figure (3);
    clf
    hp1 = pcolor(fields(1,:),freq_total,squeeze(log(imchiy)));
    set(hp1, 'edgeColor','none')
    colorbar
    title('Susceptibility in y direction')
    
% Color plot the susceptibilities
    hfig2 = figure (4);
    clf
    hp2 = pcolor(fields(1,:),freq_total,squeeze(log(imchiz)));
    set(hp2, 'edgeColor','none')
    colorbar
    title('Susceptibility in z direction')
    
%     Color plot to compare the susceptibilities of different directions
    hfig2 = figure (5);
    clf
    hp3 = pcolor(fields(1,:),freq_total,squeeze(log(imchiy)));
    hold on
    set(hp3, 'edgeColor','none')
    hp4 = pcolor(fields(1,:),freq_total,squeeze(log(imchiz)));
    set(hp4, 'edgeColor','none')
    colorbar
    
    title('Comparison of Susceptibilities')
    
% %         saveplotsT(hfig1,'im_chi')
% 
%     hfig1 = figure (6);
%     clf
%     set(hfig1,'position',[50 100 600 400])
%     hp2 = pcolor(fields(1,:),freq_total,squeeze(rechi1));
%     set(hp2, 'edgeColor','none')
%     colorbar
end