function [Hcf]=cf(B,J,Jx,Jy,Jz)

%----- Form the J operators
if nargin<3
    n=2*J+1;
    Jz=diag(J:-1:-J);
    Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
    Jm=Jp';
else
    n=size(Jx,1);
    Jp=Jx+1i*Jy;
    Jm=Jx-1i*Jy;
end


%----- Build Vcf from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X*eye(n);
O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(n);
O44c=(Jp^4+Jm^4)/2;
% O44s=-0.5i*(Jp^4-Jm^4);
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(n);

O64c=0.25*((11*Jz^2-X*eye(n)-38*eye(n))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(n)-38*eye(n)));
O64s=-1i*0.25*((11*Jz^2-X*eye(n)-38*eye(n))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(n)-38*eye(n)));

% Hcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O44s+B(5)*O60+B(6)*O64c+B(7)*O64s;
Hcf=B(1)*O20+B(2)*O40+B(3)*O44c+B(4)*O60+B(5)*O64c+B(6)*O64s;

end