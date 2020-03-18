function Ham=H3Level(hvec,Jmf)
%H0 bestimmen

% Initiate J operators
% J=8;
% Jz=diag(J:-1:-J);
% Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
% Jm=Jp';
% Jx=(Jp+Jm)/2;
% Jy=(Jp-Jm)/2i;
% 
% Ham=cf();
% [v,e]=eig(Ham);
% e=real(diag(e));
% [e,n]=sort(e);
% v=v(:,n);
% 
% Jxp=v(:,1:3)'*Jx*v(:,1:3);
% Jyp=v(:,1:3)'*Jy*v(:,1:3);
% Jzp=v(:,1:3)'*Jz*v(:,1:3);
% 
% mz=sqrt(abs(det(Jzp(1:2,1:2))));
% my=sqrt(Jyp(1:2,3)'*Jyp(1:2,3));
% mx=sqrt(Jxp(1:2,3)'*Jxp(1:2,3));
% Delta=e(3)-e(1);

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

Hcf=Delta*(eye(3)-Jz^2);

Jx=mx*Jx;
Jy=my*Jy;
Jz=mz*Jz;

Hzeeman=(-1.25*0.05788)*(hvec(1)*Jx+hvec(2)*Jy+hvec(3)*Jz);

%Jensen(2000): Jex=-0.436mueV 
Jex=-0.436;
Jex=0;
d=diag([3.912 3.912 6.821])/1000+4*Jex*eye(3)/1000;

h_dipol=d*[Jmf(1) Jmf(2) Jmf(3)]';
Ham=Hcf+Hzeeman-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;

% Include nuclear spin through hyperfine coupling?
if 1==1
  nJ=2*1+1;
  nI=2*7/2+1;
  % Initiate I operators
  Iz=diag(7/2:-1:-7/2);
  Ip=diag(sqrt((7/2-[5/2:-1:-7/2]).*(7/2+1+[5/2:-1:-7/2])),1);
  Im=Ip';
  Ix=(Ip+Im)/2;
  Iy=(Ip-Im)/2i;

  hHcfz=zeros(nJ*nI);
  hJx=zeros(nJ*nI);
  hJy=zeros(nJ*nI);
  hJz=zeros(nJ*nI);
  hIJ=zeros(nJ*nI);  
  for n=0:nI-1
    hHcfz(n*nJ+(1:nJ),n*nJ+(1:nJ))=Ham;
    hJx(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jx;
    hJy(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jy;
    hJz(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jz;
    for m=0:nI-1
      hIJ(n*nJ+(1:nJ),m*nJ+(1:nJ))=Jx*Ix(n+1,m+1)+Jy*Iy(n+1,m+1)+Jz*Iz(n+1,m+1);
    end
  end
  hIJ=0.039/11.6*hIJ;
  clear Ham
  Ham=hHcfz+hIJ;
  Jx=hJx;
  Jy=hJy;
  Jz=hJz;
end



