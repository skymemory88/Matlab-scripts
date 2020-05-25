function linear_response_theory
cd('W:\MF_calc\variousA\A097')
load('0.150.mat') % loads variables "fields", "temp", "E" and "V" 
% which are eigenstates and eigenvalues calculated in the mean-field model 
% as a function of transverse field and temperature

for zzz = 1:length(freq_total(1,:)) %calculate susceptibility for all frequencies
    freq = freq_total (zzz);
    
J=8;
I=3.5;
gLande_Ho=1.25;

%Initiate J operators
Jz=diag(J:-1:-J);
Jzh=kron(Jz,eye(2*I+1));
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jph=kron(Jp,eye(2*I+1));
Jmh=kron(Jm,eye(2*I+1));
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;
%tensor product of cristal field to include nuclear moments
%Initiate I operators
Iz=diag(I:-1:-I);
Izh=kron(eye(2*J+1),Iz);
Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im=Ip';
Iph=kron(eye(2*J+1),Ip);
Imh=kron(eye(2*J+1),Im);
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

ghztomeV = 1/241.8;
omega = freq*ghztomeV;     % define frequency sweep range (meV)
gama = 0.0001; % define lifetime (meV)

for l = 1:length(temp(1,:)) % calculate susceptibility for all temperatures
t = temp(1,l);

for k = 1:length(fields(1,:)) % calculate susceptibility for all fields
v = squeeze ( V(k,l,:,:) );
e = squeeze ( E(k,l,:) );
field = fields(1,k);
    
N = length(e);
chi_t = zeros(1,N^2);
ll = 1;
zz = zeros(1,N);
beta = 1/(t/11.6);
z=sum(exp(-beta*e));
zz=exp(-beta*e)/z;
[n,np]=meshgrid(zz,zz);
NN=n-np;
[ee,eep]=meshgrid(e,e);
EE1=1./(ee-eep-omega);
EE = eep-ee-omega;
gamma = ones(size(EE))*gama;
G = gamma ./ (EE.^2 + gamma.^2);
G1 = EE ./ (EE.^2 + gamma.^2);

ELEf = 1.250 * 0.05788;
NUCf = 4.732 * 3.1519e-5;
JxhT = Jxh * ELEf;
IxhT = Ixh * NUCf;
JyhT = Jyh * ELEf;
IyhT = Iyh * NUCf;
JzhT = Jzh * ELEf;
IzhT = Izh * NUCf;
tittt = 'S(Jyy+Iyy)';
tt  = v'  * (JyhT+IyhT) * v; 
chi_t  = (tt) .* (tt.') .* NN .* G;
chi_t1 = (tt) .* (tt.') .* NN .* G1;
sss=sum(sum(chi_t));
sss1=sum(sum(chi_t1));
imchi  (k) =  real(sss)   ;
rechi1 (k) =  real(sss1)  ;
end
        hfig1 = figure (1);
        clf
        set(hfig1,'position',[50 100 600 400])
        h1=plot (fields(1,:), imchi ,'r','LineWidth',2); 
        set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        set(gca,'fontsize',15)
        xlim([0 9]);
        xlabel('Magnetic field (T)','FontSize',15)
        ylabel('\chi'''' (arb. u.)','FontSize',15)

        hfig2 = figure (2);
        clf
        set(hfig2,'position',[680 100 600 400])
        h2=plot (fields(1,:), rechi1 ,'r','LineWidth',2); 
        set(gca,'XTick',[0,1,2,3,4,5,6,7,8,9]);
        set(gca,'fontsize',15)
        xlim([0 9]);
        xlabel('Magnetic field (T)','FontSize',15)
        ylabel('\chi'' (arb. u.)','FontSize',15)
end
end
end