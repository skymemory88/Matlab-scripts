function para_diag
clear all

Hcf=cf_Ho;
I=3.5;
J=8;
Jz=diag(J:-1:-J);
Jzh=kron(Jz,eye(2*I+1));
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jph=kron(Jp,eye(2*I+1));
Jmh=kron(Jm,eye(2*I+1));
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;
%tensor product of cristal field to include nuclear moments
Hcfh=kron(Hcf,eye(2*I+1));
%Initiate I operators
Iz=diag(I:-1:-I);
Izh=kron(eye(2*J+1),Iz);
Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im=Ip';
Iph=kron(eye(2*J+1),Ip);
Imh=kron(eye(2*J+1),Im);
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

A=0.003361;
% A_Ho_Henrik=0.003361;
% A_Ho_1985=1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000;
% A_Ho_1985=0.0066;
gLande=1.25;

h_dipol = [0.000 0.000 0.029];
Hmf = -h_dipol(1)*Jxh -h_dipol(2)*Jyh -h_dipol(3)*Jzh ;

Hhyp = A*(Ixh*Jxh+Iyh*Jyh+Izh*Jzh);
% Hhyp=A*(Izh*Jzh);
% Hsi = Hcfh + Hhyp + Hmf;
% Hsi=Hcfh;
Hsi = Hcfh + Hhyp + Hmf;


b=0:0.01:9;
B=[zeros(size(b));zeros(size(b));b];

h=1;
for h=1:size(B,2)
    Hzeeman=(-gLande*0.05788)*(B(1,h)*Jxh+B(2,h)*Jyh+B(3,h)*Jzh);
    Ham=Hsi+Hzeeman;
    [v,d]=eig(Ham);
%     e(h,:)=diag(d)-min(diag(d));
    e(h,:)=diag(d);
    
        V(h,1,:,:) = v;
        E(h,1,:) = diag(d);


    iz(h,:)=diag(v'*Izh*v);
    ix(h,:)=diag(v'*Ixh*v);
    iy(h,:)=diag(v'*Iyh*v);
    jz(h,:)=diag(v'*Jzh*v);
    jx(h,:)=diag(v'*Jxh*v);
    jy(h,:)=diag(v'*Jyh*v);
    
%     iz(:,:)=(v'*Izh*v);
%     ix(:,:)=(v'*Ixh*v);
%     iy(:,:)=(v'*Iyh*v);
%     jz(:,:)=(v'*Jzh*v);
%     jx(:,:)=(v'*Jxh*v);
%     jy(:,:)=(v'*Jyh*v);
end
hplanck = 4.13566733e-15; % eV*s  is Planck's constant
e = e / 1e3 / hplanck / 1e9 + 5100; % GHz
    
hfig4=figure(1);
plot(B(3,:),jx(:,1:16),'x')
% ylabel('Jx','FontSize',20)
% xlabel('Hc [T]','FontSize',20)
% saveplots(hfig4,'Jx')
figure(2)
plot(B(3,:),jy(:,1:16),'x')
figure(3)
hfig6=plot(B(3,:),jz(:,1:16),'x');
ylabel('Jz','FontSize',20)
xlabel('Hc [T]','FontSize',20)
% saveplots(3,'Jz')
% 
% 
hfig3=figure(4);
plot(B(3,:),ix(:,1:16),'x')
% ylabel('Ix','FontSize',20)
% xlabel('Hc [T]','FontSize',20)
% saveplots(hfig3,'Ix')
figure(5)
plot(B(3,:),iy(:,1:16),'x')
figure(6)
hfig5=plot(B(3,:),iz(:,1:16),'x');
ylabel('Iz','FontSize',20)
xlabel('Hc [T]','FontSize',20)
% saveplots(6,'Iz')

for i=1:size(B,2)
 for jj=1:15 % NMR
%     for jj=1:8 % EPR
    energies(i,jj)=e(i,jj+1)-e(i,jj); % NMR, delta_m_nuclear=1, delta_m_electron=0
%     energies(jj)=E(17-jj)-E(jj); % EPR, delta_m_nuclear=0, delta_m_electron=1
 end
end

hfig2 = figure(7);
clf
ttitle='para diag - Hc - Energy levels';
plot(B(3,:),e(:,1:16));
set(hfig2,'position',[50 50 950 950])
xlabel('Hc [T]','FontSize',20)
ylabel('E levels [GHz]','FontSize',20)
title(ttitle,'FontSize',20)
% saveplots(hfig2,ttitle)

hfig1 = figure(8);
clf
ttitle='para diag - Hc - Difference between E levels';
plot(B(3,:),energies(:,:));
set(hfig1,'position',[50 50 950 950])
xlabel('Hc [T]','FontSize',20)
ylabel('Difference between E levels [GHz]','FontSize',20)
title(ttitle,'FontSize',20)
ylim([0 5 ]);
% saveplots(hfig1,ttitle)


    fff = b;
    ttt = 0.2;
    eee = E;
    vvv = V;
        cd('W:\MF_calc\data')
        save('test','ttt','fff','eee','vvv','-v7.3')

end


function saveplots(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
curdir = cd;
cd('C:\Users\ikovacev\Documents\Projects\Plot_Matlab_Data\Plots')

saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg','pdf'});

disp(['Figure ' figname ' saved to '])
disp(cd)
cd(curdir)
end
