function chi0=chi0_w(hvec,t,momente,omega,epsilon,Jex,hyper)

if nargin<6
    Jex=0;hyper=0;
end
% Initiate J operators
J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

if hyper~=0
    I=7/2;
    Iid=eye(2*I+1);
    Iz=diag(I:-1:-I);
    Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
    Im=Ip';
    Ix=(Ip+Im)/2;
    Iy=(Ip-Im)/2i;
    hJx=kron(Jx,Iid);
    hJy=kron(Jy,Iid);
    hJz=kron(Jz,Iid);
else
    hJx=Jx;
    hJy=Jy;
    hJz=Jz;
end


d=(6/5)^2*0.05368*(dipole_direct([0 0 0],50)+exchange([0,0,0],Jex));

chi0=zeros(3,3,size(momente,1),length(omega));


for ionn=1:size(momente,1) %ionen
   h_dipol=[0 0 0];
   for ionm=1:4
        h_dipol=h_dipol+momente(ionm,:)*d(:,:,ionm,ionn)';
   end
  Ham=REHAM(hvec,h_dipol,hyper);
  
  [v,e]=eig(Ham);
  e=real(diag(e));
  e=e-min(e);
  [e,n]=sort(e);
  v=v(:,n);
  
  % calculate matric elements
  mJx=v'*hJx*v;
  mJy=v'*hJy*v;
  mJz=v'*hJz*v;
  % And corresponding weight factors:
  if t>0;
     ne=exp(-e*11.6/t)/sum(exp(-e*11.6/t));
  else
     ne=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0  
  end
  erw=[sum(diag(mJx).*ne) sum(diag(mJy).*ne) sum(diag(mJz).*ne)];
  %mJx=mJx-erw(1)*eye(size(mJx));
  %mJy=mJy-erw(2)*eye(size(mJy));
  %mJz=mJz-erw(3)*eye(size(mJz));
  
for n=1:length(omega)
w=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)))./ ...
  (e*ones(1,size(v,1))-ones(size(v,1),1)*e'-omega(n)-i*epsilon);

chi0(:,:,ionn,n)=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
      sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
      sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];
end
end
return


