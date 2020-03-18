function d=dipole(q,N,exchange)

if mod(q,1)==[0 0 0]
  warning('For q=0 or q=tau, we use the direct summation')
  d=dipole_direct(q,N);
  return
end

% Lattice parameters in aangstrom
a=[5.175 0 0
   0 5.175 0
   0 0 10.75];
% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));

% Convergence factor R~2/a
%R=2/vol^(1/3);
% Gave funne results for tau_ab~=0
R=1;

% Positions of the moments within the unit cell
tau=[0 0 0
     0 1/2 1/4
     1/2 1/2 1/2
     1/2 0 3/4];
% Convert tau to a
tau=tau*a;
%R=2/min(min(tau(tau~=0)));


b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];

q=q*b;
qq=sqrt(sum(q.*q));

%kroneckerdelta
delta=[1 0 0;0 1 0;0 0 1];


hkl=zeros((2*N+1)^3-1,3);
%hkl=[];
n=1;
for h=-N:N
  for k=-N:N
    for l=-N:N
      if h~=0 | k~=0 | l~=0
          
        %hkl=[hkl;h k l];
        hkl(n,:)=[h k l];
        n=n+1;  
      end
    end
  end
end


ql=hkl*b;

qql=[q(1)+ql(:,1) q(2)+ql(:,2) q(3)+ql(:,3)];
qqll=sum(qql.^2,2);
sumq=zeros(3,3,size(tau,1),size(tau,1));
for nt=1:size(tau,1)
for mt=1:nt
  qltau=-(tau(nt,:)-tau(mt,:))*ql';
  exp_qltau=exp(-qqll'/4/R^2-i*qltau);
  for n=1:3
  for m=1:3
    sumq(n,m,nt,mt)=sumq(n,m,nt,mt)...
            +exp_qltau*(qql(:,n).*qql(:,m)./qqll);
  end
  end
end
end

rl=hkl*a;
rll=sum(rl.^2,2);

errfc=erfc(R*sqrt(rll));
exp_qrl=exp(i*q*rl');
exp_rll=2*R*exp(-R^2*rll)/sqrt(pi)./rll;
rll15=rll.^1.5;
rll25=rll.^2.5;
sumr=zeros(3);
for n=1:3
for m=1:3
sumr(n,m)=sumr(n,m)+exp_qrl ...
          *(exp_rll ...
          .*((3+2*R^2*rll).*rl(:,n).*rl(:,m)./rll-delta(n,m)) ...
          +(3*rl(:,n).*rl(:,m)./rll25-delta(n,m)./rll15) ...
          .*errfc);
end
end

d=zeros(3,3,size(tau,1),size(tau,1));
for nt=1:size(tau,1)
for mt=1:nt
  d(:,:,nt,mt)=-4*pi/vol*(q'*q)/qq^2*exp(-qq^2/4/R^2) ...
       -4*pi/vol*sumq(:,:,nt,mt) ...
       +sumr ...
       +(4*R^3/3/sqrt(pi))*delta*(nt==mt);
  d(:,:,mt,nt)=conj(d(:,:,nt,mt));
end
end

% Sometimes the coupling is given in units of the unit cell volume:
%d=d*vol;

if nargin==3
    if exchange==0;
        return
    end
end

%Exchange via direkter Summation
clear hkl
Jex=-0.436*0.001/((5/4)^2*0.05368);
N=2;
hkl=zeros((2*N+1)^3,3);
%hkl=[];
n=1;
for h=-N:N
  for k=-N:N
    for l=-N:N
%      hkl=[hkl;h k l];
       hkl(n,:)=[h k l];
       n=n+1;
     end
  end
end


for nt=1:size(tau,1)
for mt=1:nt
  r=hkl*a;
  r(:,1)=r(:,1)-tau(nt,1)+tau(mt,1);
  r(:,2)=r(:,2)-tau(nt,2)+tau(mt,2);
  r(:,3)=r(:,3)-tau(nt,3)+tau(mt,3);
  rr=sum(r.^2,2);
  r(find(rr<0.01),:)=[];
  rr(find(rr<0.01))=[];  
  %rr15=rr.^1.5;
  rr15=rr.*sqrt(rr);
  rr25=rr.*rr15;
  exp_qr=exp(-i*q*r');
  for n=1:3
  for m=1:3
    d(n,m,nt,mt)=d(n,m,nt,mt)+exp_qr*(rr<14)*Jex*delta(n,m);%exchange
  end
  end
  d(:,:,mt,nt)=conj(d(:,:,nt,mt));
end
end

% Sometimes the coupling is given in units of the unit cell volume:
%d=d*vol;

return
