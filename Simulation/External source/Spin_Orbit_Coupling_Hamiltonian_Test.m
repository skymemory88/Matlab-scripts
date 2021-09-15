%% Generate Hso for LiErF4 / LiHoF4;
ion = LiErF4;

S = ion.S;
L = ion.L;

Xz = @(x)diag(x:-1:(-x));
Xp = @(x)diag(sqrt((x-[x-1:-1:-x]).*(x+1+[x-1:-1:-x])),1);

%----- Form the S operators
Sz=Xz(S);
Sp=Xp(S);
Sm=Sp';
Sx=(Sp+Sm)/2;
Sy=(Sp-Sm)/2i;

%----- Form the L operators
Lz=Xz(L);
Lp=Xp(L);
Lm=Lp';
Lx=(Lp+Lm)/2;
Ly=(Lp-Lm)/2i;

%% Gives the equivalent of the J-matrix in LS space

LSz = kron(Lz,eye(2*S+1));
LSx = kron(Lx,eye(2*S+1));
LSy = kron(Ly,eye(2*S+1));
LSp = LSx+1i*LSy;
LSm = LSx-1i*LSy;

SLz = kron(eye(2*L+1),Sz);
SLx = kron(eye(2*L+1),Sx);
SLy = kron(eye(2*L+1),Sy);
SLp = SLx+1i*SLy;
SLm = SLx-1i*SLy;

% Constants from Christensen 1979, effective hamiltonian from Karayianis 1970, and Karayanis and Farrar 1970
LdotS = LSz*SLz+LSx*SLx+LSy*SLy;
Hso = ion.lambda(1)*LdotS+ion.lambda(2)*LdotS^2+ion.lambda(3)*LdotS^3;
[e_so,v_so] = diag_ham(Hso,0);
[e_so_unique,j2p1_so] = unique(round(e_so,4));
j2p1_so = [j2p1_so;length(e_so)+1];

%% Crystal Fields in LS-space
X = @(J)(J*(J+1)); %J = L+S, so it works still. 

%----- Form the Operators
    
Z = LSz+SLz;
P = LSp+SLp;
M = LSm+SLm;

O20_LS=3*Z^2-X(J)*eye((2*L+1)*(2*S+1));
O40_LS=35*Z^4-(30*X(J)-25)*Z^2+(3*X(J)^2-6*X(J))*eye((2*L+1)*(2*S+1));
O44_LS=(P^4+(P')^4)/2;
O60_LS=231*Z^6-(315*X(J)-735)*Z^4+(105*X(J)^2-525*X(J)+294)*Z^2+(-5*X(J)^3+40*X(J)^2-60*X(J))*eye((2*L+1)*(2*S+1));
O64c_LS=0.25*((11*Z^2-X(J)*eye((2*L+1)*(2*S+1))-38*eye((2*L+1)*(2*S+1)))*(P^4+(P')^4)+(P^4+(P')^4)*(11*Z^2-X(J)*eye((2*L+1)*(2*S+1))-38*eye((2*L+1)*(2*S+1))));
O64s_LS=-1i*0.25*((11*Z^2-X(J)*eye((2*L+1)*(2*S+1))-38*eye((2*L+1)*(2*S+1)))*(P^4-(P')^4)+(P^4-(P')^4)*(11*Z^2-X(J)*eye((2*L+1)*(2*S+1))-38*eye((2*L+1)*(2*S+1))));

Hcf_LS=B_LS(1)*O20_LS+B_LS(2)*O40_LS+B_LS(3)*O44_LS+B_LS(4)*O60_LS+B_LS(5)*O64c_LS+B_LS(6)*O64s_LS;
[e_cf_so] = diag_ham(Hcf_LS+Hso,0);
%         line([0 1],e_cf_J(1)*ones(2,1),'Color','k','LineStyle','--')
%     for i = 2:length(e_cf_J)
%         if round(e_cf_J(i-1),4)==round(e_cf_J(i),4)
%         line([0 1],e_cf_J(i)*ones(2,1),'Color','b','LineStyle','--')
%         else
%         line([0 1],e_cf_J(i)*ones(2,1),'Color','k','LineStyle','--')
%         end
%     end
%     set(gca,'TickLabelInterpreter','latex','FontSize',14)
%     xticklabels({})
%     ylabel('Energy (meV)','Interpreter','latex')

%% Adding in the HF term

I = ion.I;
Iz=Xz(I);
Ip=Xp(I);
Im=Ip';
Ix=(Ip+Im)/2;
Iy=(Ip-Im)/2i;

Hsoh = kron(Hso,eye(2*I+1));
Hcfh = kron(Hcf_LS,eye(2*I+1));
Zh=kron(Z,eye(2*I+1));
Ph=kron(P,eye(2*I+1));
Mh=kron(M,eye(2*I+1));
Xh=(Ph+Mh)/2;
Yh=(Ph-Mh)/2i;

Izh=kron(eye((2*L+1)*(2*S+1)),Iz);
Iph=kron(eye((2*L+1)*(2*S+1)),Ip);
Imh=kron(eye((2*L+1)*(2*S+1)),Im);
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

ELEf = ion.gLande*0.05788;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173

%Calculate Hamiltonian
hvec = [0 0 0];
Hzeeman = -ELEf*(hvec(1)*Xh+hvec(2)*Yh+hvec(3)*Zh);
H_hyper = ion.A*(Ixh*Xh + Iyh*Yh + Izh*Zh);

Ham = Hcfh+Hzeeman+H_hyper+Hsoh;
[e_total] = diag_ham(Ham,1);

%% Spectre output
eig_Spectre = unique([0 0 2.1 2.1 3.4 3.4 7.7 7.7 38.4 38.4 43.9 43.9 48.2 48.2 51.5 51.5,...
815.9 815.9 816.7 816.7 822.5 822.5 838 838 840.9 840.9 845.4 845.4 847.2 847.2,...
1273.7 1273.7 1276.3 1276.3 1284.9 1284.9 1288 1288 1288.7 1288.7 1290.5 1290.5,...
1534.3 1534.3 1557.4 1557.4 1563.6 1563.6 1566.7 1566.7 1581.2 1581.2,...
1907.8 1907.8 1911.5 1911.5 1914.1 1914.1 1924.9 1924.9 1932.8 1932.8,...
2298.6 2298.6 2307.7 2307.7,...
2391.7 2391.7 2394.8 2394.8 2400.1 2400.1 2410.9 2410.9 2411.5 2411.5 2413.3 2413.3,... 
2562 2562 2562.1 2562.1 2576.3 2576.3 2578.2 2578.2,... 
2773.8 2773.8 2775.4 2775.4 2781.5 2781.5 2818.2 2818.2 2823.4 2823.4]);

eig_Christensen = [0 16 21 50 232 273 305 334,...
6785 6789 6819 6907 6938 6956 6975,...
10243 10256 10296 10312 10322 10329,...
11927 12054 12102 12123 12218,...
15334 15354 15365 15429 15479,...
18319 18380 18912 18938 18972 19008 19004 19025,...
20657 20663 20754 20768,...
22310 22318 22356,...
22553 22578,...
24019 24142 24163 24188 24273]*0.12398;

figure;
    subplot(2,2,1)
    plot(h,squeeze(e_so_cf_hf_zee(:,1,:)),'ok')
    xlabel('H$\parallel$a (T)','Interpreter','latex')
    ylabel('Energy (meV)','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    ylim([0 50])
    subplot(2,2,2)
    plot(h,squeeze(e_so_cf_hf_zee(:,1,:)),'ok')
    xlabel('H$\parallel$a (T)','Interpreter','latex')
    ylabel('Energy (meV)','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    ylim([800 880])
    subplot(2,2,3)
    plot(h,squeeze(e_so_cf_hf_zee(:,1,:)),'ok')
    xlabel('H$\parallel$a (T)','Interpreter','latex')
    ylabel('Energy (meV)','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    ylim([1250 1300])
    subplot(2,2,4)
    plot(h,squeeze(e_so_cf_hf_zee(:,1,:)),'ok')
    xlabel('H$\parallel$a (T)','Interpreter','latex')
    ylabel('Energy (meV)','Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    ylim([1550 1600])

e0 = e_so;
e1 = eig_Spectre;
e2 = squeeze(e_cf_so);
e3 = diag_ham(cf(ion.J,ion.B),0);
    figure;
    subplot(4,2,1:2:8)
        hold on
        for i = 1:length(e1)
            line([0 1],e1(i)*ones(2,1),'Color','k')
        end
        for i = 1:length(e2)
            line([0 1],e2(i)*ones(2,1),'Color','b')
        end
        for i = 1:length(e3)
            line([0 1],e3(i)*ones(2,1),'Color','r','LineStyle','--')
        end
        for i = 1:length(e0)
            line([0 1],e0(i)*ones(2,1),'Color','g','LineStyle','--','LineWidth',2)
        end
        set(gca,'TickLabelInterpreter','latex','FontSize',14)
        xticklabels({})
        ylabel('Energy (meV)','Interpreter','latex')
    subplot(4,2,8)
        hold on
        for i = 1:length(e1)
            line([0 1],e1(i)*ones(2,1),'Color','k')
        end
        for i = 1:length(e2)
            line([0 1],e2(i)*ones(2,1),'Color','b')
        end
        for i = 1:length(e3)
            line([0 1],e3(i)*ones(2,1),'Color','r','LineStyle','--')
        end
        for i = 1:length(e0)
            line([0 1],e0(i)*ones(2,1),'Color','g','LineStyle','--')
        end
        set(gca,'TickLabelInterpreter','latex','FontSize',14)
        xticklabels({})
        ylabel('Energy (meV)','Interpreter','latex')
        ylim([0 60])
    subplot(4,2,6)
        hold on
        for i = 1:length(e1)
            line([0 1],e1(i)*ones(2,1),'Color','k')
        end
        for i = 1:length(e2)
            line([0 1],e2(i)*ones(2,1),'Color','b')
        end
        for i = 1:length(e3)
            line([0 1],e3(i)*ones(2,1),'Color','r','LineStyle','--')
        end
        for i = 1:length(e0)
            line([0 1],e0(i)*ones(2,1),'Color','g','LineStyle','--')
        end
        set(gca,'TickLabelInterpreter','latex','FontSize',14)
        xticklabels({})
        ylabel('Energy (meV)','Interpreter','latex')
        ylim([800 850])
    subplot(4,2,4)
        hold on
        for i = 1:length(e1)
            line([0 1],e1(i)*ones(2,1),'Color','k')
        end
        for i = 1:length(e2)
            line([0 1],e2(i)*ones(2,1),'Color','b')
        end
        for i = 1:length(e3)
            line([0 1],e3(i)*ones(2,1),'Color','r','LineStyle','--')
        end
        for i = 1:length(e0)
            line([0 1],e0(i)*ones(2,1),'Color','g','LineStyle','--')
        end
        set(gca,'TickLabelInterpreter','latex','FontSize',14)
        xticklabels({})
        ylabel('Energy (meV)','Interpreter','latex')
        ylim([1250 1300])
    subplot(4,2,2)
        hold on
        for i = 1:length(e1)
            line([0 1],e1(i)*ones(2,1),'Color','k')
        end
        for i = 1:length(e2)
            line([0 1],e2(i)*ones(2,1),'Color','b')
        end
        for i = 1:length(e3)
            line([0 1],e3(i)*ones(2,1),'Color','r','LineStyle','--')
        end
        for i = 1:length(e0)
            line([0 1],e0(i)*ones(2,1),'Color','g','LineStyle','--')
        end
        set(gca,'TickLabelInterpreter','latex','FontSize',14)
        xticklabels({})
        ylabel('Energy (meV)','Interpreter','latex')
        ylim([1500 1600])        
