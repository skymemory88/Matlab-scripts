function lin_resp_th_works_omegascan_A097
close all
% temp=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.5 3 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95 4];
temp=[0.15];

% temp=[0.1];
% freq = [3.436,   0.4674,   1.682, 3.330, 5.038,   0.85, 2.55, 4.2, 5.9,   3.925,   4.45];
% freq = [1.682 3.436 3.924 4.449 5.604 0.000001];
freq = [1.682];
for j = 1:length(temp(1,:))
lname=[num2str(temp(j),'%3.3f'),'.mat']; 
ik_read_absorption_f(lname,freq)
end

% temp=0.20;
% freq = 3.42;
% % freq = 3.925;
% for j = 1:length(temp(1,:))
% lname=[num2str(temp(j)),'.mat'];
% ik_read_absorption_f(lname,freq)
% end

end

function ik_read_absorption_f(lname,freq_total)
clearvars -except lname freq_total 
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption\data')
cd('W:\MF_calc\dataA097\data')
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption\data\A-hyp\1')
% lname='0.02.mat';
% save('fieldscan.mat','temp','fields','ion','history','E','V')
% save(tit,'ttt','fff','eee','vvv','-v7.3')
load(lname)
% load('0.200-5degrees.mat')

for zzz = 1:length(freq_total(1,:))
    freq = freq_total (zzz);
    

J=8;
I=3.5;
gLande_Ho=1.25;

% With hyperfine coupling
%Initiate J operators
Jz=diag(J:-1:-J);
Jzh=kron(Jz,eye(2*I+1));
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
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
% epsilon=0.0004;                % define lineshape width (meV)
% sigma=0.0003; % define lineshape width
gama = 0.0001;

fields=fff;
temp=ttt;
E=eee;
V=vvv;
% clear field 

for l = 1:length(temp(1,:))
t = temp(1,l);

for k = 1:length(fields(1,:))

    v = squeeze ( V(k,l,:,:) );
    e = squeeze ( E(k,l,:) );
    field = fields(1,k);
    
N = length(e);
chi_t = zeros(1,N^2);
ll = 1;
zz = zeros(1,N);

beta = 1/(t/11.6);
z=sum(exp(-beta*e));
zz=exp(-beta*e)/z;
[n,np]=meshgrid(zz,zz);
NN=n-np;
[ee,eep]=meshgrid(e,e);
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

tittt = 'S(Jyy+Iyy)';
tt  = v'  * (JyhT+IyhT) * v; 
% tt  = v'  * (cos(pi/180*10)*JyhT + sin(pi/180*10)*JzhT) * v;

chi_t  = (tt) .* (tt.') .* NN .* G;
chi_t1 = (tt) .* (tt.') .* NN .* G1;

% chi_el_tt = diag(tt) ;
% chi_el_ss = sum (chi_el_tt .* chi_el_tt .* zz ) ;
% chi_el_av = sum (chi_el_tt .* zz) ;
% chi_el = beta * (chi_el_ss - chi_el_av * chi_el_av) ;


sss=sum(sum(chi_t));
sss1=sum(sum(chi_t1));

imchi  (k) =  real(sss)   ;
% imchi  (k) =  real(sss) +  real(chi_el) *  omega*gama /(omega^2 + gama^2)  ;
% imchi  (k) =  imag(sss );
rechi1 (k) =  real(sss1)  ;
% rechi1 (k) =  real(sss1) + real(chi_el) *  gama^2 /(omega^2 + gama^2)  ;
% rechi1 (k) =  imag(sss1);
% S11   (k) =  abs (1i*real(sss) + real(sss1));
% elastic   (k) =  real (chi_el);

end
        hfig1 = figure (1);
        clf
%         hold on
        set(hfig1,'position',[50 100 600 400])
%         if freq == 3.436
        h1=plot (fields(1,:), imchi ,'r','LineWidth',2); 
%         else
%         h1=plot (fields(1,:), imchi ,'b','LineWidth',2); 
%         end
        set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        set(gca,'fontsize',15)
        xlim([0 9]);
%         ylim([0 1.1]);
        xlabel('Magnetic field (T)','FontSize',15)
        ylabel('\chi'''' (arb. u.)','FontSize',15)
%         tit=num2str(t,omega,'Temperature = %3.3f K, omega = %d '); 
        tit=[tittt,', \chi'''', T = ',num2str(t,'%3.3f'),' K, \omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, \Gamma = ',num2str(gama,'%3.5f'),' meV'];
        tit1=[tittt,', chi'''', T = ',num2str(t,'%3.3f'),' K, omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, Gamma = ',num2str(gama,'%3.5f'),' meV'];
        title(tit,'FontSize',10)
        grid off
        box on
%         sname=num2str(t,'Imaginary %3.3f K');
%         saveplotsT(hfig1,tit1)


        hfig2 = figure (2);
        clf
%         hold on
        set(hfig2,'position',[680 100 600 400])
%         if freq == 3.436
        h2=plot (fields(1,:), rechi1 ,'r','LineWidth',2); 
%         else
%         h2=plot (fields(1,:), rechi1 ,'b','LineWidth',2); 
%         end
        set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        set(gca,'fontsize',15)
        xlim([0 9]);
%         ylim([0 1.1]);
        xlabel('Magnetic field (T)','FontSize',15)
        ylabel('\chi'' (arb. u.)','FontSize',15)
%         tit=num2str(t,omega,'Temperature = %3.3f K, omega = %d '); 
        tit=[tittt,', \chi'', T = ',num2str(t,'%3.3f'),' K, \omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, \Gamma = ',num2str(gama,'%3.5f'),' meV'];
        tit1=[tittt,', chi'', T = ',num2str(t,'%3.3f'),' K, omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, Gamma = ',num2str(gama,'%3.5f'),' meV'];
        title(tit,'FontSize',10)
        grid off
        box on
%         sname=num2str(t,'Imaginary %3.3f K');
%         saveplotsT(hfig2,tit1)

%         hfig3 = figure (3);
%         clf
% %         hold on
%         set(hfig3,'position',[1300 100 600 400])
%         if freq == 3.436
%         h2=plot (fields(1,:), elastic ,'r','LineWidth',2); 
%         else
%         h2=plot (fields(1,:), elastic ,'b','LineWidth',2); 
%         end
%         set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
%         set(gca,'fontsize',15)
%         xlim([0 9]);
% %         ylim([0 1.1]);
%         xlabel('Magnetic field (T)','FontSize',15)
%         ylabel('elastic (arb. u.)','FontSize',15)
% %         tit=num2str(t,omega,'Temperature = %3.3f K, omega = %d '); 
%         tit=[tittt,', elastic, T = ',num2str(t,'%3.3f'),' K, \omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, \Gamma = ',num2str(gama,'%3.5f'),' meV'];
%         tit1=[tittt,', elastic, T = ',num2str(t,'%3.3f'),' K, omega = ',num2str(omega/ghztomeV,'%3.3f'),' GHz, Gamma = ',num2str(gama,'%3.5f'),' meV'];
%         title(tit,'FontSize',10)
%         grid off
%         box on
% %         sname=num2str(t,'Imaginary %3.3f K');
%         saveplotsT(hfig3,tit1)
end
end

%         hfig1 = figure (3);
%         clf
%         hp1 = pcolor(fields(1,:),freq_total,log(imchi)');
%         set(hp1, 'edgeColor','none')
%         colorbar
%         saveplotsT(hfig1,'im_chi')
%         
%         
%         hfig1 = figure (4);
%         clf
%         set(hfig1,'position',[50 100 600 400])
%         hp2 = pcolor(fields(1,:),freq_total,rechi1');
%         set(hp2, 'edgeColor','none')
%         colorbar
%         saveplotsT(hfig1,'re_chi')

end

function saveplotsT(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
cd('D:\Projects\LiHoF4_Network_Analyzer\plots\new_analysis\A097\lin_resp_th')

saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg'});

disp(['Figure ' figname ' saved to '])
disp(cd)
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption')
end