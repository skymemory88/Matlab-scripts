function [y, name, pnames, pin]=cfpfit(x, p, flag)

if nargin==2
x=x(:)';
faktor=p(1);
B=p(2:8);
    
az=49;
ber=1:49;

env1=x(ber);
ber=ber+az;
qxx1=x(ber);
ber=ber+az;
qzz1=x(ber);
ber=ber+az;
qxz1=x(ber);
ber=ber+az;
q1=x(ber);
ber=ber+az;
int1=x(ber);
ber=ber+az;
env2=x(ber);
ber=ber+az;
qxx2=x(ber);
ber=ber+az;
qzz2=x(ber);
ber=ber+az;
qxz2=x(ber);
ber=ber+az;
q2=x(ber);
ber=ber+az;
int2=x(ber);
ber=ber+az;
env3=x(ber);
ber=ber+az;
qxx3=x(ber);
ber=ber+az;
qzz3=x(ber);
ber=ber+az;
qxz3=x(ber);
ber=ber+az;
q3=x(ber);
ber=ber+az;
int3=x(ber);
    
    
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

Jx12=v(:,1:2)'*Jx*v(:,3:4);
Jz12=v(:,1:2)'*Jz*v(:,3:4);  
Jx13=v(:,1:2)'*Jx*v(:,5:6);
Jz13=v(:,1:2)'*Jz*v(:,5:6);  
Jx14=v(:,1:2)'*Jx*v(:,7:8);
Jz14=v(:,1:2)'*Jz*v(:,7:8);  

%formfactor:

ff1=formfactor(q1).^2;
ff2=formfactor(q2).^2;
ff3=formfactor(q3).^2;

I1=abs(faktor*ff1.*( qxx1*sum(sum(Jx12'.*Jx12))+qzz1*sum(sum(Jz12'.*Jz12))+2*qxz1*sum(sum(Jx12'.*Jz12))  ))
I2=abs(faktor*ff2.*( qxx2*sum(sum(Jx13'.*Jx13))+qzz2*sum(sum(Jz13'.*Jz13))+2*qxz2*sum(sum(Jx13'.*Jz13))  ))
I3=abs(faktor*ff3.*( qxx3*sum(sum(Jx14'.*Jx14))+qzz3*sum(sum(Jz14'.*Jz14))+2*qxz3*sum(sum(Jx14'.*Jz14))  ))

en1=(e(3)-e(1))*ones(1,az);
en2=(e(5)-e(1))*ones(1,az);
en3=(e(7)-e(1))*ones(1,az);

y=[en1,qxx1,qzz1,qxz1,q1, I1,...
   en2,qxx2,qzz2,qxz2,q2, I2,...
   en3,qxx3,qzz3,qxz3,q3, I3]';







else
   y=[];
   name='cfpfit';
   pnames=str2mat('Amplitude','Centre','Width');
	
end
return

function f=formfactor(q)
f=0*q+1;


