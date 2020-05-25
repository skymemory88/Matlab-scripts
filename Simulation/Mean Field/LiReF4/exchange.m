function [d]=exchange(q,Jex)

% This function performs a brute force summation of
% the q-dependent exchange coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% nr is the number of unit cells that should be summed
% in each direction.

% Parameters for LiHoF4
% Unit vectors in aangstroms
a=[5.162 0 0
   0 5.162 0
   0 0 10.70];
% Positions of the moments within the unit cell
tau=[0 0 0
     0 1/2 1/4
     1/2 1/2 1/2
     1/2 0 3/4];
%Parameters for Ho2Ti2O7
%a=[10.12 0 0
%   0 10.12 0
%   0 0 10.12];
%tau=[0   0   0
%     0   1/2 1/2
%     1/2 0 1/2
%     1/2 1/2 0
%     0   1/4 1/4
%     0   3/4 3/4
%     1/2 1/4 3/4
%     1/2 3/4 1/4
%     1/4 1/4 0
%     1/4 3/4 1/2
%     3/4 1/4 1/2
%     3/4 3/4 0
%     1/4 0   1/4
%     1/4 1/2 3/4
%     3/4 0   3/4
%     3/4 1/2 1/4];

% Convert tau to a
tau=tau*a;

% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];
% Convert q from Miller indicies to reciprocal aangstroms
q=q*b;
% Length of q
qq=sqrt(sum(q.*q));

% Kronecker delta in x,y,and z
delta=[1 0 0;0 1 0;0 0 1];

% For N moments in the unit cell, there will be N 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij=J_ji), so we just calculate J_1j. 
% The result is a (3x3xN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.

N=1;

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


d=zeros(3,3,size(tau,1),size(tau,1));
for nt=1:size(tau,1)
for mt=1:nt
  r=hkl*a;
  r(:,1)=r(:,1)-tau(nt,1)+tau(mt,1);
  r(:,2)=r(:,2)-tau(nt,2)+tau(mt,2);
  r(:,3)=r(:,3)-tau(nt,3)+tau(mt,3);
  rr=sum(r.^2,2);
  r(find(rr<0.01),:)=[];
  rr(find(rr<0.01))=[];  
  r(find(rr>(N*5.175)^2),:)=[];%%%%%% clear Spins outside of the sphere
  rr(find(rr>(N*5.175)^2))=[]; %%%%%%
  exp_qr=exp(-1i*q*r');
  for n=1:3
  for m=1:3
    d(n,m,nt,mt)=exp_qr*(rr<14)*Jex*delta(n,m);%exchange
  end
  end
  %  d(:,:,nt,mt)=d(:,:,nt,mt)+(4*pi/3)*0.01389*eye(3)/4; %Lorentz
  d(:,:,mt,nt)=conj(d(:,:,nt,mt));
end
end

% Sometimes the coupling is given in units of the unit cell volume:
%d=d*vol;


return
