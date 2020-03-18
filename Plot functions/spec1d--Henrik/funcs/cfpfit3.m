function [y, name, pnames, pin]=cfpfit3(x, p, flag)

if nargin==2

global qglobal  

qxx=qglobal(:,1);
qxz=qglobal(:,2);
qzz=qglobal(:,3);
ff2=qglobal(:,4);
ni=qglobal(:,5);
nf=qglobal(:,6);
temp=qglobal(:,7);


b0=[-0.678 -6.82 -0.1332 -0.0243];
%B=[67.7 -0.678 -6.82 -0.00008 -0.1332 -0.0243]/1000;
%xb = fsolve(@(xb)ediff1(xb,[p(1:2)]),b0);
%xb=b0;


%B=[p(1) xb(1) xb(2) p(2) xb(3) xb(4) ]/1000;
B=p/1000;
B*1000
xpt=[1:size(qglobal,1)]';
ypt=0*xpt;


[mxx,myy,mzz,mxz,ei]=intnm(B);

for n=1:length(xpt)
    z=exp(-11.6*ei(ni(n))/temp(n))/sum(exp(-11.6*ei/temp(n)));
    ypt(n)=ff2(n)*z*(qxx(n)*mxx(ni(n),nf(n))+qxz(n)*mxz(ni(n),nf(n))+qzz(n)*mxz(ni(n),nf(n))+ myy(ni(n),nf(n)));
%
end

y=p(7)*ypt;

else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Cfpfit';
	pnames=str2mat('B20','B64c','B64s','Skalierung');
	if flag==1, pin=[0 0 0 0]; else pin = p; end
   
end

return 

function [mxx,myy,mzz,mxz,ei]=intnm(B)

%----- Form the J operators

J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;
dimension=2*J+1;

%----- Build Vcf from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X*eye(16);
O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(16);
O44c=(Jp^4+Jm^4)/2;
%O44s=-1i*(Jp^4-Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);
O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

Vcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O60+B(5)*O64c+B(6)*O64s;%+B(4)*O44s

%----- calculate Eigenvalues

[v,e]=eig(Vcf);
e=real(diag(e));
[e,n]=sort(e);
v=v(:,n);
e=e-min(e);

MJx=v'*Jx*v;
MJz=v'*Jz*v;
MJy=v'*Jy*v;

jxjx=MJx.*conj(MJx);
jxjz=MJx.*conj(MJz)+MJz.*conj(MJx);
jzjz=MJz.*conj(MJz);
jyjy=MJy.*conj(MJy);

%Uebergang Dublett ni nach Dublett mf

for ni=1:dimension/2
    ei(ni)=e(2*ni-1);
    for mf=1:dimension/2
        mber=(2*mf-1):(2*mf);
        nber=(2*ni-1):(2*ni);
        mxx(ni,mf)=sum(sum(jxjx(nber,mber)));
        mzz(ni,mf)=sum(sum(jzjz(nber,mber)));
        mxz(ni,mf)=sum(sum(jxjz(nber,mber)));
        myy(ni,mf)=sum(sum(jyjy(nber,mber)));
    end
end

return

%Energielevels justieren:

function y=ediff1(x,Brt)
global eglobal
%B=[ Brt(1)  x(1) x(2) Brt(2) x(3) x(4)]/1000;
B=x/1000;
%----- Form the J operators

J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

%----- Build Vcf from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X*eye(16);
O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(16);
O44c=(Jp^4+Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);
O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

B=[67.7 -0.678 -6.82 -0.00008 -0.1332 -0.0243]/1000;
Vcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O60+B(5)*O64c+B(6)*O64s;

e=real(eig(Vcf));

y=e([3 5 7])-e(1)-eglobal';




%%%%%%%%%%%%%%%Fragmente:

%T=fliplr(eye(2*J+1));

%O64s*T+T*O64s

%----- calculate Eigenvalues



%pv1=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
%pv2=[0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
%pv3=[0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
%pv4=[0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
%Ts=[eye(4) zeros(4) zeros(4) zeros(4);
%    zeros(4) eye(4) eye(4) eye(4)]


%T=eye(2*J+1);
%T=T(:,[2:2*J+1 1]);

%T*O64c-conj(O64c)*T
%P=pv1'*pv1; 
 
%red1=pv1*Vcf*pv1'
%red1=pv2*Vcf*pv2'
%red1=pv3*Vcf*pv3'
%red1=pv4*Vcf*pv4'

%P*Vcf-Vcf*P

%sop1(:,:,1)=pv1*O20*pv1';
%sop1(:,:,2)=pv1*O40*pv1';
%sop1(:,:,3)=pv1*O44c*pv1';
%sop1(:,:,4)=pv1*O60*pv1';
%sop1(:,:,5)=pv1*O64c*pv1';
%sop1(:,:,6)=pv1*O64s*pv1';

%sop2(:,:,1)=pv2*O20*pv2';
%sop2(:,:,2)=pv2*O40*pv2';
%sop2(:,:,3)=pv2*O44c*pv2';
%sop2(:,:,4)=pv2*O60*pv2';
%sop2(:,:,5)=pv2*O64c*pv2';
%sop2(:,:,6)=pv2*O64s*pv2';

%sop3(:,:,1)=pv3*O20*pv3';
%sop3(:,:,2)=pv3*O40*pv3';
%sop3(:,:,3)=pv3*O44c*pv3';
%sop3(:,:,4)=pv3*O60*pv3';
%sop3(:,:,5)=pv3*O64c*pv3';
%sop3(:,:,6)=pv3*O64s*pv3';



%h1=zeros(4);
%h2=zeros(4);
%h3=zeros(4);
%for n=1:6
    %h1=h1+B(n)*sop1(:,:,n);
    %h2=h2+B(n)*sop2(:,:,n);
 %   h3=h3+B(n)*sop3(:,:,n);
%end

%eig(h1)
%eig(h2)-eig(h3)


%es=sort([eig(h1)', eig(h2)']);

%y=es(2:4)-min(es)-[2.2 3.5 7];


