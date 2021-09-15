function [Ham,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] =calc_ham(ion,hvec,h_dipol,ishf)

J = ion.J;
I = ion.I;

Jz=diag(J:-1:-J);
Jzh=kron(Jz,eye(2*I+1));
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jph=kron(Jp,eye(2*I+1));
Jmh=kron(Jm,eye(2*I+1));
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;

%tensor product of crystal field to include nuclear moments
Hcf = cf(ion.J,ion.B);
Hcfh=kron(Hcf,eye(2*I+1));

%Initiate I operators
Iz=diag(I:-1:-I);
Izh=kron(eye(2*J+1),Iz);
Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im=Ip';
Iph=kron(eye(2*J+1),Ip);
Imh=kron(eye(2*J+1),Im);
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

ELEf = ion.gLande*0.05788;     % Lande factor * Bohr magneton (meV T^-1)
NUCf = 4.173 * 3.15245e-5;   % Nuclear Lande factor, mu/mu_N = 4.173

%Calculate Hamiltonian
Hzeeman = -ELEf*(hvec(1)*Jxh+hvec(2)*Jyh+hvec(3)*Jzh);
HzeemanI = -NUCf*(hvec(1)*Ixh+hvec(2)*Iyh+hvec(3)*Izh);
H_dipol = -h_dipol(1)*Jxh-h_dipol(2)*Jyh-h_dipol(3)*Jzh;
H_hyper = ion.A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh);

Ham = Hcfh + Hzeeman + HzeemanI + H_dipol + H_hyper;


end
