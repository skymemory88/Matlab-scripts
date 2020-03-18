function [y, name, pnames, pin]=strfit(x, p, flag)


if nargin==2
    global hklglobal
    y=0*x;
    for n=1:size(hklglobal,1)
        [Fs,e,lorentz]=strfm(hklglobal(n,:),p(2:length(p)));
        Fsp=Fs-(e*Fs)*e';
        y(n)= abs(p(1))*real(Fsp'*Fsp)/lorentz;
    end
else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Power law (y=0 x>x0)';
	pnames=str2mat('skalierung','X.offset','Exponent','Background');
	if flag==1, pin=[0 0 0 0]; else pin = p; end
   if flag==2
        %-------------- Get starting guesses 
        % Get starting guess for x offset
        mf_msg('Click on x offset estimate (y=0 x>x0)');	    
        [xoff dummy]=ginput(1);

        % Get starting guess for background 
        mf_msg('Click on background estimate');	    
        [dummy bg]=ginput(1);

        % Get two points on the 'curve'
        mf_msg('Click on two points on curve');	    
        [xp yp]=ginput(2);			   

        % Work out starting guess for the exponent
        power=log((yp(1)-bg)/(yp(2)-bg))/log((xp(1)-xoff)/(xp(2)-xoff));

        % Work out starting guess for the amplitude
        amp=abs((yp(1)-bg)/((xp(1)-xoff)^power));
    
        pin=[amp xoff power bg];
    end
end

return

function [F,e,lorentz]=strfm(hkl,mp)

a=[5.182 0 0
   0 5.182 0
   0 0 10.7];

tau=[    0    0.2500    0.6250
    0.5000    0.2500    0.8750
    0.5000    0.7500    0.1250
    0.0000    0.7500    0.3750];


% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];

% Convert q from Miller indicies to reciprocal aangstroms
q=hkl*b;
% Length of q
qq=sqrt(sum(q.*q,2));
lorentz=1.18*qq/(4*pi);%lorentzfaktor

e=q./qq;

ff2=formfactorq(qq);

%spin=[1 1 1 1;
%      0 0 0 0;
%      0 0 0 0];
%spin(:,1)=mp(2)*spin(:,1);

global modellglobal

spin=zeros(3,4);
for nmp=1:length(mp)
    spin=spin+mp(nmp)*modellglobal(:,:,nmp);
end
    
F=zeros(3,1);
for nt=1:4
    F=F+ff2*spin(:,nt)*exp(-2i*pi*hkl*tau(nt,:)');
end

return

function f=formfactorq(q)
L=6;
S=3/2;
J=15/2;

A=0.0389;   
a=5.3118;   
B=0.2598;   
b=8.1732;   
C=0.6784;   
c=2.0828;  
D=0.0222;


s=sqrt(q.*q)/(4*pi);
j0=A*exp(-a*s.^2)+B*exp(-b*s.^2)+C*exp(-c*s.^2)+D;
j1=(s.^2).*j0;

f=j0+j1*(J*(J+1)-S*(S+1)+L*(L+1))/(3*J*(J+1)+S*(S+1)-L*(L+1));

return