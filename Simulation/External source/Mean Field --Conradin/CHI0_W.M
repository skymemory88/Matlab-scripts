function chi0=chi0_w(hvec,t,Jmf,omega,epsilon)

mx=3.3517;
my=mx;
mz=5.5111;
Delta=0.9345;

% Initiate J operators
J=1;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

Jx=mx*Jx;
Jy=my*Jy;
Jz=mz*Jz;

% Include nuclear spin through hyperfine coupling?
if 1==1
  nJ=2*1+1;
  nI=2*7/2+1;
  
  hJx=zeros(nJ*nI);
  hJy=zeros(nJ*nI);
  hJz=zeros(nJ*nI);
  hIJ=zeros(nJ*nI);  
  for n=0:nI-1
    hJx(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jx;
    hJy(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jy;
    hJz(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jz;
  end
  Jx=hJx;
  Jy=hJy;
  Jz=hJz;
end

Ham=H3Level(hvec,Jmf);
[v,e]=eig(Ham);
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);

% calculate matric elements
mJx=v'*Jx*v;
mJy=v'*Jy*v;
mJz=v'*Jz*v;

% And corresponding weight factors:
if t>0;
  ne=exp(-e*11.6/t)/sum(exp(-e*11.6/t));
else
  ne=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0  
end

for n=1:length(omega)
w=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)))./ ...
  (e*ones(1,size(v,1))-ones(size(v,1),1)*e'-omega(n)-i*epsilon);

chi0(:,:,n)=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
      sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
      sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];
end
return


