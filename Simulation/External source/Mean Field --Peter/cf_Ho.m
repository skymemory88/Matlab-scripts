function [Hcf]=cf_Ho()


%----- Form the J operators

J=8;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

%----- Build Vcf from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X*eye(17);
O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(17);
O44=(Jp^4+Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(17);

O64c=0.25*((11*Jz^2-X*eye(17)-38*eye(17))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(17)-38*eye(17)));
O64s=-i*0.25*((11*Jz^2-X*eye(17)-38*eye(17))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(17)-38*eye(17)));

% Parameters in meV according to JJ's note
% B=[-60 0.35 3.6 0.0004 0.07 0.0098]/1000; % in meV
B=[-60 0.35 3.6 0.0004 0.07 0.006]/1000; % in meV

%B=[-64.8 .4258 4.5344 .0001003 .08553 .01688]/1000; % in meV

Hcf=B(1)*O20+B(2)*O40+B(3)*O44+B(4)*O60+B(5)*O64c+B(6)*O64s;
