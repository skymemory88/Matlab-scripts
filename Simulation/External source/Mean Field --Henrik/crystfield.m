function [split,vg,ve,Hcfz,v,e]=crystfield(h)

% Below is the crystal field operator as
% reconstructed from the table VII in
% P.E. Hansen et al. PRB vol. 12 p. 5315 (1975)

% Energies are in Kelvin -> meV (1 meV=11.6 K)
%ecf=[535 492 398 92 74 546 546 497 497 124 124 0 0 575 511 39 9]/11.6;
ecf=[535 492 398 92 74 546 497 124 0 546 497 124 0 575 511 39 9]/11.6;

cs1=[exp(i*1.7*pi/180) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 exp(-i*8.2*pi/180) 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 exp(i*8.2*pi/180) 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 exp(-i*1.7*pi/180)];

cw1=[0.397 -.437 .551 -.437 .397;...
 .53 -.468 0 .468 -.530;...
-.386 .158 .807 .158 -.386;...
.468 .530 0 -.530 -.468;...
.440 .533 .212 .533 .440];

cs2a=[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
      0 0 0 0 0 exp(-i*9.2*pi/180) 0 0 0 0 0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0 0 0 exp(-i*88.1*pi/180) 0 0 0 0 0 0 0;...
      0 0 0 0 0 0 0 0 0 0 0 0 0 exp(-i*92.5*pi/180) 0 0 0];
cs2b=fliplr(conj(cs2a));

cw2=[.562 -.717 -.371 .183;...
     .281 -.305 .788 -.456;...
     .070 .023 .488 .870;...
     .775 .627 -.060 -.045];

cs3=[0 0 exp(i*30.4*pi/180) 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 exp(i*22.7*pi/180) 0 0 0 0 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 exp(-i*22.7*pi/180) 0 0 0 0 0 0;...
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 exp(-i*30.4*pi/180) 0 0];

cw3=[.374 -.600 -.600 .374;...
     -.423 .566 -.566 .423;...
     .600 .374 .374 .600;...
     -.566 -.423 .423 .566];

cf=[cw1*cs1;cw2*cs2a;cw2*cs2b;cw3*cs3];

Hcf=cf'*diag(ecf)*cf;

%----- Form the J operators

J=8;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

%----- Build Vcf from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X;
O40=35*Jz^4-(30*X-25)*Jz^2+3*X^2-6*X;
O44=(Jp^4+Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2-5*X^3+40*X^2-60*X;
Or64=(11*Jz^2-X-38)*Jp^4+Jp^4*(11*Jz^2-X-38);
Or6m4=(11*Jz^2-X-38)*Jm^4+Jm^4*(11*Jz^2-X-38);
O64c=sqrt( 1/2)*(Or64+Or6m4);
O64s=sqrt(-1/2)*(Or6m4-Or64);

O64c=0.25*((11*Jz^2-X-38)*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X-38));

% Parameters in meV according to JJ's note
B=[-64.8 .4258 4.5344 .0001003 .08553 .01688]/1000; % in meV

%Hcf=B(1)*O20+B(2)*O40+B(3)*O44+B(4)*O60+B(5)*O64c; %+B(6)*O64s;

%----- Zeeman term from transverse field:
% 1 Bohr magneton = 0.6717  Kelvin / Tesla
%                 = 0.05788 meV / Tesla
% Holmium Lande factor is 5/4.

%Hzeeman=(-5/4*0.05788*h)*Jx;
% The transverse field g-factor is 0.74,
% But I do not see why I should use this
% like David did.
if length(h)==1
  Hzeeman=(-0.74*0.05788*h)*Jx;
else
  Hzeeman=(-0.74*0.05788)*(h(1)*Jx+h(2)*Jy+h(3)*Jz);
end
  
Hcfz=Hcf+Hzeeman;

% Diagonalize to get energy levels and eigenvectors
[v,e]=eig(Hcfz);
[e,n]=sort(real(diag(e)));
v=v(:,n);
% e contains the eigenvalues.
% the rows of v contains the eigenvectors.

split=e(2)-e(1);
vg=v(:,1);
ve=v(:,2);

% To plot the levels and splitting as a function of
% applied field, comment out the `return; above.

ecf=[];
h=0:0.1:8;
for n=1:length(h)
%  Hzeeman=-5/4*0.05788*h(n)*Jx;
  Hzeeman=-0.74*0.05788*h(n)*Jx;
  [ecf(:,n),nn]=sort(real(eig(Hcf+Hzeeman)));
  [evec,eval]=eig(Hcf+Hzeeman);
  jzg(n)=evec(:,nn(1))'*Jz*evec(:,nn(1));
  jze(n)=evec(:,nn(2))'*Jz*evec(:,nn(2));
  jzeg(n)=evec(:,nn(2))'*Jz*evec(:,nn(1));
  jzge(n)=evec(:,nn(1))'*Jz*evec(:,nn(2));
  jxg(n)=evec(:,nn(1))'*Jx*evec(:,nn(1));
  jxe(n)=evec(:,nn(2))'*Jx*evec(:,nn(2));
  jxeg(n)=evec(:,nn(2))'*Jx*evec(:,nn(1));
  jxge(n)=evec(:,nn(1))'*Jx*evec(:,nn(2));
  jyg(n)=evec(:,nn(1))'*Jy*evec(:,nn(1));
  jye(n)=evec(:,nn(2))'*Jy*evec(:,nn(2));
  jyeg(n)=evec(:,nn(2))'*Jy*evec(:,nn(1));
  jyge(n)=evec(:,nn(1))'*Jy*evec(:,nn(2));
end

jzg=abs(jzg);
jze=abs(jze);
jzeg=abs(jzeg);
jzge=abs(jzge);

%figure(1)
%plot(h,jzg,h,jze,h,jzeg,h,jzge)
%figure(2)
%plot(h,jxg,h,jxe,h,jxeg,'x-',h,jxge)
%figure(3)
%plot(h,jyg,h,jye,h,jyeg,h,jyge)

return

if 1==1
split=ecf;

plot(h,split)
axis([0 8 -2 2.5])
set(gcf,'paperposition',[.5 .5 5 3.75])
set(gca,'position',[.2 .15 .7 .7],'fontname','times','fontsize',20)
title('Lowest crystal field levels in LiHoF_4')
xlabel('Transverse field [T]')
ylabel('Energy [meV]')
print -depsc crystfield.eps
end

gmub=0.74*0.05788;
clf
ha=plot(h,split(2,:)-split(1,:),h,0.5*gmub*h.*jxg,'r--')
set(ha,'linewidth',2)
axis([0 8 0 1])
set(gcf,'paperposition',[.5 .5 5 3.75])
set(gca,'position',[.2 .15 .7 .7],'fontname','times','fontsize',20)
xlabel('Transverse field [T]')
ylabel('\Gamma [meV]')

set(gca,'position',[.15 .15 .8 .8])
axes('position',[.17 .53 .4 .4])
ha=plot(h,split)
set(ha,'linewidth',2)
axis([0 8 -1 1.5])
set(gca,'fontname','times','fontsize',16)
ylabel('Energy [meV]')
set(gca,'yaxislocation','right')
xlabel('Transverse field [T]')
print -depsc -loose lihof4_ising1.eps

%break

clf
jz=real([jzg(1) jzge(2:end)]);
ha=plot(h,jz,h,jxg,'r--')
set(ha,'linewidth',2)
%axis([0 8 -2 2.5])
set(gcf,'paperposition',[.5 .5 5 3.75])
set(gca,'position',[.15 .15 .8 .8],'fontname','times','fontsize',20)
xlabel('Transverse field [T]')
ylabel('M, A')
print -depsc -loose lihof4_ising2.eps

%return