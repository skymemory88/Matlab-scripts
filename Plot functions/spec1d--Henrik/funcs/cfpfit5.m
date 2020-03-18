function [y, name, pnames, pin]=cfpfit4(x, p, flag)

if nargin==2

global qglobal izst fzst temperaturen skalierung 

J=15/2;
B=p(1:7)/1000;
y=0*x;

[matrel,eif,en]=intnm(B,J);

for ny=1:length(izst)
    besetzung=exp(-11.6*en(izst(ny))/temperaturen(ny))/sum(exp(-11.6*en/temperaturen(ny)));
    y(ny)=real(p(skalierung(ny))*besetzung*sum(sum( qglobal(:,:,ny).*squeeze(matrel(izst(ny),fzst(ny),:,:)) )));
    if izst(ny)==1&fzst(ny)==3
        besetzung=exp(-11.6*en(3)/temperaturen(ny))/sum(exp(-11.6*en/temperaturen(ny)));
        y(ny)=y(ny)+real(p(skalierung(ny))*besetzung*sum(sum( qglobal(:,:,ny).*squeeze(matrel(3,4,:,:)) )));
    end    
end

y(ny+1:ny+3)=en(2:4);
y(ny+4:ny+6)=en([5,6,8])/10;
%y(ny+5)=[36/10]';
%y(ny+4)=32/10;
%y(ny+6)=43.5/10;
en
%Suszeptibilitaet;

[chixx2, chizz2]=suszept(1000,B);
faktor=10000*0.058/(6.0228e23*9.274e-21);%molar
y(end)=faktor./chizz2-faktor./chixx2;




% 
% temperatur=2;
% 
% [mxx,myy,mzz,mxz,en]=intnm(B);
% 
% y=0*x;
% =en(2:4)-en(1);
% 
% 
% y(4)= p(7)*z(1)*( qsumme(1,1)*mxx(1,2)+qsumme(1,2)*mxz(1,2)+qsumme(1,3)*mzz(1,2)+qsumme(1,4)*myy(1,2) );
% y(5)= p(7)*z(1)*( qsumme(2,1)*mxx(1,3)+qsumme(2,2)*mxz(1,3)+qsumme(2,3)*mzz(1,3)+qsumme(2,4)*myy(1,3) );
% y(6)= p(7)*z(1)*( qsumme(3,1)*mxx(1,4)+qsumme(3,2)*mxz(1,4)+qsumme(3,3)*mzz(1,4)+qsumme(3,4)*myy(1,4) );
% temperatur=20;
% z=exp(-11.6*en/temperatur)/sum(exp(-11.6*en/temperatur));
% y(7)= p(8)*z(1)*( qsumme(4,1)*mxx(1,2)+qsumme(4,2)*mxz(1,2)+qsumme(4,3)*mzz(1,2)+qsumme(4,4)*myy(1,2) );
% y(8)= p(8)*z(1)*( qsumme(5,1)*mxx(1,3)+qsumme(5,2)*mxz(1,3)+qsumme(5,3)*mzz(1,3)+qsumme(5,4)*myy(1,3) );
% y(9)= p(8)*z(1)*( qsumme(6,1)*mxx(1,4)+qsumme(6,2)*mxz(1,4)+qsumme(6,3)*mzz(1,4)+qsumme(6,4)*myy(1,4) );
% y(10)=p(8)*z(2)*( qsumme(7,1)*mxx(2,3)+qsumme(7,2)*mxz(2,3)+qsumme(7,3)*mzz(2,3)+qsumme(7,4)*myy(2,4) );
% y(11)=p(8)*z(2)*( qsumme(8,1)*mxx(2,4)+qsumme(8,2)*mxz(2,4)+qsumme(8,3)*mzz(2,4)+qsumme(8,4)*myy(2,3) );
% 
% 
% %[chixx1, chizz1]=suszept(900,B);
% [chixx2, chizz2]=suszept(1000,B);
% 
% faktor=10000*0.058/(6.0228e23*9.274e-21);%molar
% 
% %y(7)=faktor./chizz1-faktor./chixx1;
% y(12)=faktor./chizz2-faktor./chixx2;

else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Cfpfit';
	pnames=str2mat('B20','B64c','B64s','Skalierung');
	if flag==1, pin=[0 0 0 0]; else pin = p; end
   
end

return 

function [matrel,eif,en]=intnm(B,J)

global O20 O40 O44s O44c O60 O64s O64c Jx Jy Jz

%----- Form the J operators

J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;
dimension=2*J+1;

%----- Build Vcf from the stevens operators:

% X=J*(J+1);
% O20=3*Jz^2-X*eye(16);
% O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(16);
% O44c=(Jp^4+Jm^4)/2;
% O44s=-1i*(Jp^4-Jm^4)/2;
% O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);
% O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
% O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

Vcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O60+B(5)*O64c+B(6)*O64s+B(7)*O44s;

%----- calculate Eigenvalues

[v,e]=eig(Vcf);
e=real(diag(e));
[e,n]=sort(e);
v=v(:,n);
e=e-min(e);
en=e(1:2:end);
eif=en*ones(size(en'))-ones(size(en))*en';
%Uebergang Dublett ni nach Dublett mf

MJ=zeros(2*J+1,2*J+1,3);    %ni, nf,alpha
MJ(:,:,1)=v'*Jx*v; 
MJ(:,:,2)=v'*Jy*v;
MJ(:,:,3)=v'*Jz*v;

matrel=zeros(J+0.5,J+0.5,3,3);

Projektor=  [  1     0     0     0     0     0     0     0
               1     0     0     0     0     0     0     0
               0     1     0     0     0     0     0     0
               0     1     0     0     0     0     0     0
               0     0     1     0     0     0     0     0
               0     0     1     0     0     0     0     0
               0     0     0     1     0     0     0     0
               0     0     0     1     0     0     0     0
               0     0     0     0     1     0     0     0
               0     0     0     0     1     0     0     0
               0     0     0     0     0     1     0     0
               0     0     0     0     0     1     0     0
               0     0     0     0     0     0     1     0
               0     0     0     0     0     0     1     0
               0     0     0     0     0     0     0     1
               0     0     0     0     0     0     0     1 ];

for alpha=1:3
    for beta=1:3
        matrel(:,:,alpha,beta)=Projektor'*(MJ(:,:,alpha).*conj( MJ(:,:,beta)  )*Projektor);
    end
end

return

function [chixx, chizz]=suszept(t,B)

%----- Form the J operators

global O20 O40 O44s O44c O60 O64s O64c Jx Jy Jz

J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;
dimension=2*J+1;

%----- Build Vcf from the stevens operators:
% 
% X=J*(J+1);
% O20=3*Jz^2-X*eye(16);
% O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(16);
% O44c=(Jp^4+Jm^4)/2;
% O44s=-1i*(Jp^4-Jm^4)/2;
% O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);
% O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
% O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

Vcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O60+B(5)*O64c+B(6)*O64s+B(7)*O44s;

%----- calculate Eigenvalues

[v,e]=eig(Vcf);
e=real(diag(e));
[e,n]=sort(e);
v=v(:,n);
e=e-min(e);

g=1.2;
muB=0.05788;

e=e-min(e);
z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));



vjzv=abs(v'*Jz*v).^2;
vjxv=abs(v'*Jx*v).^2;

nenner=e*ones(1,2*J+1)-ones(2*J+1,1)*e';
nnull=(abs(nenner)<0.1);
nenner=(ones(2*J+1)-nnull).*nenner+(t/11.6)*nnull;

pij=(ones(2*J+1)-nnull).*(-z*ones(1,2*J+1)+ones(2*J+1,1)*z')+nnull*diag(z);

chixx=(g*muB)^2*sum(sum(vjxv.*pij./nenner));
chizz=(g*muB)^2*sum(sum(vjzv.*pij./nenner));
return
