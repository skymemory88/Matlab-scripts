function energies_lowest3_Ha_hyperfine
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

Hhyp=A*(Ixh*Jxh+Iyh*Jyh+Izh*Jzh);
% Hhyp=A*(Izh*Jzh);
Hsi=Hcfh+Hhyp;
% Hsi=Hcfh;

b=0:0.01:9;
B=[b;zeros(size(b));zeros(size(b))];

for h=1:size(B,2)
    Hzeeman=(-gLande*0.05788)*(B(1,h)*Jxh+B(2,h)*Jyh+B(3,h)*Jzh);
    Ham=Hsi+Hzeeman;
    [v,d]=eig(Ham);
%     e(h,:)=diag(d)-min(diag(d));
if h==1
zero=min(diag(d));
end
e(h,:)=diag(d)-zero;
%     e(h,:)=diag(d);
    iz(h,:)=diag(v'*Izh*v);
    ix(h,:)=diag(v'*Ixh*v);
    iy(h,:)=diag(v'*Iyh*v);
    jz(h,:)=diag(v'*Jzh*v);
    jx(h,:)=diag(v'*Jxh*v);
    jy(h,:)=diag(v'*Jyh*v);
end
hplanck = 4.13566733e-15; % eV*s  is Planck's constant
e = e / 1e3 / hplanck / 1e9 ; % GHz
    
% hfig4=figure(1);
% plot(B(3,:),jx(:,1:16),'x')
% % ylabel('Jx','FontSize',20)
% % xlabel('Hc [T]','FontSize',20)
% % saveplots(hfig4,'Jx')
% figure(2)
% plot(B(3,:),jy(:,1:16),'x')
% figure(3)
% hfig6=plot(B(3,:),jz(:,1:16),'x');
% ylabel('Jz','FontSize',20)
% xlabel('Hc [T]','FontSize',20)
% saveplots(3,'Jz')
% 
% 
% hfig3=figure(4);
% plot(B(3,:),ix(:,1:16),'x')
% % ylabel('Ix','FontSize',20)
% % xlabel('Hc [T]','FontSize',20)
% % saveplots(hfig3,'Ix')
% figure(5)
% plot(B(3,:),iy(:,1:16),'x')
% figure(6)
% hfig5=plot(B(3,:),iz(:,1:16),'x');
% ylabel('Iz','FontSize',20)
% xlabel('Hc [T]','FontSize',20)
% saveplots(6,'Iz')
% 
% for i=1:size(B,2)
%  for jj=1:15 % NMR
% %     for jj=1:8 % EPR
%     energies(i,jj)=e(i,jj+1)-e(i,jj); % NMR, delta_m_nuclear=1, delta_m_electron=0
% %     energies(jj)=E(17-jj)-E(jj); % EPR, delta_m_nuclear=0, delta_m_electron=1
%  end
% end

hfig2 = figure(9);
clf
ttitle='energies_lowest3_Ha_hyperfine_zoom';
% ttitle='energies_lowest3_Ha_hyperfine';
p1=plot(B(1,:),e(:,:),'k');

set(hfig2,'position',[50 50 600 450])
% xlabel('Ha [T]','FontSize',13)
% ylabel('Frequency [GHz]','FontSize',13)
% ylabel('Relative energies [K]','FontSize',20)
% title(ttitle,'FontSize',20)
ylim([-50 50]);
xlim([0 2])
xlabel('Transverse magnetic field (T)','FontSize',13);
ylabel('Frequency (GHz)','FontSize',13)
ax=gca;
set(ax,'FontSize',13);
box on
set(p1,'linewidth',1,'Color',[35 107 142]/255,'linestyle','-')


% saveplots(hfig2,ttitle)


end


% function saveplots(hfig,figname)
% % if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
% curdir = cd;
% cd('D:\Projects\LiHoF4_Network_Analyzer\plots\thesis_figures\introduction')
% 
% saveas(figure(hfig),[figname '.fig'],'fig');
% print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
% print2eps(figname,hfig)
% [result,msg] = eps2xxx([figname '.eps'],{'jpeg','pdf'});
% 
% disp(['Figure ' figname ' saved to '])
% disp(cd)
% cd(curdir)
% end
