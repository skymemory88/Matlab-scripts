function [d]=dipole_r(hkl)

% This function performs a brute force summation of
% the q-dependent dipole coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% nr is the number of unit cells that should be summed
% in each direction.

%global 
global_exchange=0;
Jex=global_exchange;%-0.001/(6/5)^2/0.05368;

% Parameters for LiHoF4
% Unit vectors in aangstroms
a=[5.175 0 0
   0 5.175 0
   0 0 10.75];
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

% Kronecker delta in x,y,and z
delta=[1 0 0;0 1 0;0 0 1];

% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];

  r=hkl*a;
  rr=sum(r.^2,2);
  rr15=rr.*sqrt(rr);
  rr25=rr.*rr15;
  q=[0 0 0];
  exp_qr=exp(-i*q*r');
  
  for n=1:3
  for m=1:3
    d(n,m)=exp_qr*(3*r(:,n).*r(:,m)./rr25-delta(n,m)./rr15);
    d(n,m)=d(n,m)+exp_qr*(rr<14)*Jex*delta(n,m);%exchange
  end
  end
%  d(:,:,nt,mt)=d(:,:,nt,mt)+(4*pi/3)*0.01389*eye(3)/4; %Lorentz



return