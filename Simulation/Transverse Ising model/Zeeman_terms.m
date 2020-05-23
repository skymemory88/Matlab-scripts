function Zeeman_terms

mpwd = 13.9e-3;

% molar mass
molarmass = 245.43;

% number of magnetic ions/formula
Nions = 1;

% number of mols of magnetic ions
Nmols = Nions*(mpwd/molarmass);

% measuring field
Hmsr = 100;    % in Oe

LiDyF4 = importdata('G:\My Drive\Programming scripts\Matlab\0066.dc.dat',',',31);

xpwd = LiDyF4.data(:,4);
ypwd = LiDyF4.data(:,5)/Hmsr/Nmols; % in emu/mol Oe

%-----------------------file reading complete--------------------------


% Crystal field parameter coefficients:
B = [-130 0.65 7.2 0 -0.055 -0.0084]/1e3; %[meV] converted from Romanova, Mag Res Sol 8, 1 (2006)
S = 5/2;
L = 5;
J = S + L;
muB = 0.0579625;
gJ = 3/2 + (S*(S+1) - L*(L+1))/(2*J*(J+1));
k = 0.0862875;

% Define operators:
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;
n=2*J+1;

% 1. ----- Build Crystal Field Hamiltonian, Hcf, from the stevens operators:

X=J*(J+1);
O20=3*Jz^2-X*eye(n);
O40=35*Jz^4-(30*X-25)*Jz^2+(3*X^2-6*X)*eye(n);
O44=(Jp^4+Jm^4)/2;
O60=231*Jz^6-(315*X-735)*Jz^4+(105*X^2-525*X+294)*Jz^2+(-5*X^3+40*X^2-60*X)*eye(n);

O64c=0.25*((11*Jz^2-X*eye(n)-38*eye(n))*(Jp^4+Jm^4)+(Jp^4+Jm^4)*(11*Jz^2-X*eye(n)-38*eye(n)));
O64s=-1i*0.25*((11*Jz^2-X*eye(n)-38*eye(n))*(Jp^4-Jm^4)+(Jp^4-Jm^4)*(11*Jz^2-X*eye(n)-38*eye(n)));

Hcf=B(1)*O20+B(2)*O40+B(3)*O44+B(4)*O60+B(5)*O64c+B(6)*O64s;

Scal = 6.0;
H = 0.001;
Hzez = -muB*gJ*Jz*H;
Hzex = -muB*gJ*Jx*H;
Hzey = -muB*gJ*Jy*H;

Htotz = Hcf + Hzez;
Htotx = Hcf + Hzex;
Htoty = Hcf + Hzey;

Temp = xpwd';


%ajz = zeros(size(Temp));
xsiz = zeros(1,length(Temp));
xsix = zeros(1,length(Temp));
xsiy = zeros(1,length(Temp));
%xsizz = zeros(size(Temp));


for j = 1:length(Temp);
       
    Beta = 1/(k*Temp(j));
    
    [Vz, Ez] = eig(Htotz);
    shift = min(diag(Ez));
    Ez = Ez - shift.*eye(n);
    [Vx, Ex] = eig(Htotx);
    [Vy, Ey] = eig(Htoty);
    
    Zz = trace(exp(-Beta*Ez));
    Zx = trace(exp(-Beta*Ex));
    Zy = trace(exp(-Beta*Ey));
    xsiz(j) = muB*gJ*trace(real(Vz'*Jz*Vz).*exp(-Beta.*Ez)/Zz)/H;
    xsix(j) = muB*gJ*trace(real(Vx'*Jx*Vx).*exp(-Beta.*Ex)/Zx)/H;
    xsiy(j) = muB*gJ*trace(real(Vy'*Jy*Vy).*exp(-Beta.*Ey)/Zy)/H;
    
    % ajz(j) = real(diag(Vz'*Jz*Vz)'*(exp(-Beta.*diag(Ez))/Zz));
    % xsizz(j) = muB*gJ*ajz(j)/H;
    
    j = j + 1;
end

N = zeros(1,length(Temp)); % Demagnitizing factor
for j = 1:length(Temp);
    N(j) = (xsix(j)-Temp(j))/(xsix(j)*Temp(j));
end

% There is some problem with demagnetizing factor

%plot(Temp, 1./xsiz,'s', Temp,1./xsix,'o',Temp, 1./xsiy,'-')
subplot(2,2,1), 
plot(Temp, xsiz,'s', Temp, xsix,'o',Temp, xsiy,'+', xpwd, ypwd, 'x')
title('Magnetc susceptability')
subplot(2,2,2),
plot(Temp, 1./xsiz,'s', Temp, 1./xsix,'o',Temp, 1./xsiy,'-', xpwd, 1./ypwd, 'x')
title('Inverted magnetc susceptability')
subplot(2,2,3),
plot(Temp, Scal*xsiz,'s', Temp, Scal*xsix,'o',Temp, Scal*xsiy,'+', xpwd, ypwd, 'x')
title('Scaled magnetc susceptability')
subplot(2,2,4),
plot(Temp, N)
title('Demagnetizing factor')
%plot(Temp, ajz, 's')
%plot(Temp, 1./ajz, 's')
end