function [y, name, pnames, pin]=cfpfit(x, p, flag)

if nargin==2
    
anzset=length(x)/Anzahlbins/2; %Anzahlbins Golbale Variable
y=[];
B=p(1:7);
faktor=p(8);
[mxx,mxz,mzz,Eif]=intnm(B)
for n=1:anzset   
    [qxx, qxz, qzz, f22, ni, mf]=variblen(n)
    I=ff2.*(qxx*mxx(ni,mf)+qxz*mxz(ni,mf)+qzz*mxz(ni,mf));
    energ=Eif(ni,mf)*ones(1,Anzahlbins);
    y=[y,I,energ];
end

        
        
function [mxx,mxz,mzz,E]=intnm(B)

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
O44s=-1i*(Jp^4-Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);
O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

Vcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O44s+B(5)*O60+B(6)*O64c+B(7)*O64s;

%----- calculate Eigenvalues

[v,e]=eig(Vcf);
e=real(diag(e));
[e,n]=sort(e);
v=v(:,n);

MJx=v'*Jx*v;
MJz=v'*Jx*v;

jxjx=MJx.*MJx';
jxjz=MJx.*MJz'+MJz.*MJx';
jzjz=MJz.*MJz';

%Uebergang Dublett ni nach Dublett mf

for ni=1:dimension/2
    for mf=1:dimension/2
        mber=(2*mf-1):(2*mf);
        nber=(2*ni-1):(2*ni);
        mxx(ni,mf)=sum(sum(jxjx(nber,mber)));
        mzz(ni,mf)=sum(sum(jzjz(nber,mber)));
        mxz(ni,mf)=sum(sum(jxjz(nber,mber)));
    end
end

return

%formfaKtor:

function f=formfactor(q)
f=0*q+1;
return

