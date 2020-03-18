function [Hcf]=cf_Er()


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
O44=(Jp^4+Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(16);

O64c=0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));
O64s=-i*0.25*((11*Jz^2-X*eye(16)-38*eye(16))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(16)-38*eye(16)));

% Parameters in meV 

%according to Hansen
%B=[67.7 -0.678 -6.82 -0.00008 -0.1332 -0.0243]/1000; % in meV

%new fit C.Kraemer
B=[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]/1000;
%B=[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   0]/1000;
%B=[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   0]/1000;



Hcf=B(1)*O20+B(2)*O40+B(3)*O44+B(4)*O60+B(5)*O64c+B(6)*O64s;
