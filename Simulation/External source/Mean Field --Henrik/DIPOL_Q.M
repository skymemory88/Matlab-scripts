if 1==0
%----- Convergence of direct and Ewald summations respectively:
ni=1:15;
de=zeros(3,3,4,4,length(ni));
dd=zeros(3,3,4,4,length(ni));
for n=1:length(ni)
  de(:,:,:,:,n)=dipole([0.5,0,0],ni(n));
  dd(:,:,:,:,n)=dipole_direct([0.5,0,0],ni(n));
end
vol=287.8917;
de=de*vol;
dd=dd*vol;
nt=1;mt=1;
m=3;n=3;
plot(ni,squeeze(abs(de(n,m,nt,mt,:))),...
     ni,squeeze(abs(dd(n,m,nt,mt,:))),'--')
set(gca,'position',[.2 .15 .7 .7],...
        'fontname','times','fontsize',20)
axis([1 15 1.8 2.6])
set(gca,'ytick',[2 2.5])
title('Direct and Ewald summation for dipole coupling')
xlabel('Number of unit cell shells')
ylabel('Dipole coupling')
legend('Ewald''s method','Direct summation',4)
set(gcf,'paperposition',[.5 .5 6 4.5])
print -depsc dipole_convergence.eps
end % Convergence


%q=[0 0.01:0.01:0.09 0.1:0.1:2];
q=0:0.1:1;
qz=[0 0.01 0.05:0.05:0.95 0.99 1 1.01 1.05:0.05:1.95 1.99 2]'*[0 0 1];
qy=[0 0.01 0.05:0.05:0.95 0.99 1]'*[0 1 0];
qx=[0 0.01 0.05:0.05:0.95 0.99 1]'*[1 0 0];

dqx=zeros(3,3,4,4,length(qx));
dqy=zeros(3,3,4,4,length(qy));
dqz=zeros(3,3,4,4,length(qz));
for n=1:length(qz)
  dqz(:,:,:,:,n)=dipole(qz(n,:),12);
%  ddqz(:,:,:,:,n)=dipole_direct(qz(n,:),8);
end
for n=1:length(qx)
  dqx(:,:,:,:,n)=dipole(qx(n,:),12);
%  ddqx(:,:,:,:,n)=dipole_direct(qx(n,:),8);
end
for n=1:length(qx)
  dqy(:,:,:,:,n)=dipole(qy(n,:),12);
%  ddqy(:,:,:,:,n)=dipole_direct(qy(n,:),8);
end

vol=287.8917;
dqx=dqx*(5/4)^2*0.05368/vol;
dqy=dqy*(5/4)^2*0.05368/vol;
dqz=dqz*(5/4)^2*0.05368/vol;

save dipole_saves

n=1;m=1;
h=plot(qz(:,3),squeeze(real(dqz(n,m,1,1,:))),...
     qz(:,3),squeeze(real(dqz(n,m,1,2,:))),'--',...
     qz(:,3),squeeze(real(dqz(n,m,1,3,:))),'-.',...
     qz(:,3),squeeze(real(dqz(n,m,1,4,:))),':',...
     qz(:,3),sum(squeeze(real(dqz(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{xx} [meV]')
xlabel('(0,0,q) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(2.1,dqz(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(2.3,sum(dqz(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
axis([0 2 -.005 .006])
print -depsc dipole_qz_xx.eps

n=2;m=2;
h=plot(qz(:,3),squeeze(real(dqz(n,m,1,1,:))),...
     qz(:,3),squeeze(real(dqz(n,m,1,2,:))),'--',...
     qz(:,3),squeeze(real(dqz(n,m,1,3,:))),'-.',...
     qz(:,3),squeeze(real(dqz(n,m,1,4,:))),':',...
     qz(:,3),sum(squeeze(real(dqz(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{yy} [meV]')
xlabel('(0,0,q) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(2.1,dqz(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(2.3,sum(dqz(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
axis([0 2 -.005 .006])
print -depsc dipole_qz_yy.eps

n=3;m=3;
h=plot(qz(:,3),squeeze(real(dqz(n,m,1,1,:))),...
     qz(:,3),squeeze(real(dqz(n,m,1,2,:))),'--',...
     qz(:,3),squeeze(real(dqz(n,m,1,3,:))),'-.',...
     qz(:,3),squeeze(real(dqz(n,m,1,4,:))),':',...
     qz(:,3),sum(squeeze(real(dqz(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{zz} [meV]')
xlabel('(0,0,q) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(2.1,dqz(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(2.1,dqz(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(2.3,dqz(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(2.1,sum(dqz(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
axis([0 2 -.01 .006])
print -depsc dipole_qz_zz.eps


n=1;m=1;
h=plot(qx(:,1),squeeze(real(dqx(n,m,1,1,:))),...
     qx(:,1),squeeze(real(dqx(n,m,1,2,:))),'--',...
     qx(:,1),squeeze(real(dqx(n,m,1,3,:))),'-.',...
     qx(:,1),squeeze(real(dqx(n,m,1,4,:))),':',...
     qx(:,1),sum(squeeze(real(dqx(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{xx} [meV]')
xlabel('(q,0,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqx(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.15,sum(dqx(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qx_xx.eps

n=2;m=2;
h=plot(qx(:,1),squeeze(real(dqx(n,m,1,1,:))),...
     qx(:,1),squeeze(real(dqx(n,m,1,2,:))),'--',...
     qx(:,1),squeeze(real(dqx(n,m,1,3,:))),'-.',...
     qx(:,1),squeeze(real(dqx(n,m,1,4,:))),':',...
     qx(:,1),sum(squeeze(real(dqx(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{yy} [meV]')
xlabel('(q,0,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqx(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.15,dqx(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.05,sum(dqx(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qx_yy.eps

n=3;m=3;
h=plot(qx(:,1),squeeze(real(dqx(n,m,1,1,:))),...
     qx(:,1),squeeze(real(dqx(n,m,1,2,:))),'--',...
     qx(:,1),squeeze(real(dqx(n,m,1,3,:))),'-.',...
     qx(:,1),squeeze(real(dqx(n,m,1,4,:))),':',...
     qx(:,1),sum(squeeze(real(dqx(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{zz} [meV]')
xlabel('(q,0,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqx(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqx(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.15,dqx(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.05,sum(dqx(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qx_zz.eps


n=1;m=1;
h=plot(qy(:,2),squeeze(real(dqy(n,m,1,1,:))),...
     qy(:,2),squeeze(real(dqy(n,m,1,2,:))),'--',...
     qy(:,2),squeeze(real(dqy(n,m,1,3,:))),'-.',...
     qy(:,2),squeeze(real(dqy(n,m,1,4,:))),':',...
     qy(:,2),sum(squeeze(real(dqy(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{xx} [meV]')
xlabel('(0,q,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqy(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.15,dqy(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.05,sum(dqy(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qy_xx.eps

n=2;m=2;
h=plot(qy(:,2),squeeze(real(dqy(n,m,1,1,:))),...
     qy(:,2),squeeze(real(dqy(n,m,1,2,:))),'--',...
     qy(:,2),squeeze(real(dqy(n,m,1,3,:))),'-.',...
     qy(:,2),squeeze(real(dqy(n,m,1,4,:))),':',...
     qy(:,2),sum(squeeze(real(dqy(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{yy} [meV]')
xlabel('(0,q,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqy(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.15,sum(dqy(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qy_yy.eps

n=3;m=3;
h=plot(qy(:,2),squeeze(real(dqy(n,m,1,1,:))),...
     qy(:,2),squeeze(real(dqy(n,m,1,2,:))),'--',...
     qy(:,2),squeeze(real(dqy(n,m,1,3,:))),'-.',...
     qy(:,2),squeeze(real(dqy(n,m,1,4,:))),':',...
     qy(:,2),sum(squeeze(real(dqy(n,m,1,:,:)))),'-');
nicefig
set(h,'linewidth',2')
title('Dipole parameters in LiHoF_4')
ylabel('J^{zz} [meV]')
xlabel('(0,q,0) [rlu]')
set(gca,'position',[.2 .15 .6 .7])
text(1.05,dqy(n,m,1,1,end-1),'J_{11}','fontname','times','fontsize',20)
text(1.15,dqy(n,m,1,2,end-1),'J_{12}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,3,end-1),'J_{13}','fontname','times','fontsize',20)
text(1.05,dqy(n,m,1,4,end-1),'J_{14}','fontname','times','fontsize',20)
text(1.05,sum(dqy(n,m,1,:,end-1)),'sum','fontname','times','fontsize',20)
%axis([0 2 -.01 .006])
print -depsc dipole_qy_zz.eps

