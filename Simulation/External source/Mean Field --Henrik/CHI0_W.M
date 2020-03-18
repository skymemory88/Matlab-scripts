function chi0=chi0_w(h,t,omega,epsilon,lambda_dipole)

if nargin<5
  lambda_dipole=1;
end

% This function

% Get eigenvalues and vectors
if lambda_dipole==1;
  [Ham,jx,jy,jz,e,v]=meanfield(h,t,0,0,1);
else
  [split,vg,ve,Ham,v,e]=crystfield(h);
end

% Initiate J operators
J=8;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

% calculate matric elements
mJx=v'*Jx*v;
mJy=v'*Jy*v;
mJz=v'*Jz*v;

% And corresponding weight factors:
if t>0;
  ne=exp(-e*11.6/t)/sum(exp(-e*11.6/t));
else
  ne=(e==min(e));
end

w=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)))./ ...
  (e*ones(1,size(v,1))-ones(size(v,1),1)*e'-omega-i*epsilon);

chi0=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
      sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
      sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];

return


