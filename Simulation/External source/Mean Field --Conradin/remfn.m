function [jx,jy,jz]=remfn(h,t,jx,jy,jz)

% This program performs a simple meanfield calculation
% for LiHoF4

%------------Mf berechnen


% Initiate J operators
J=8;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

d=diag([3.912 3.912 6.821])/1000-4*0.436*eye(3)/1000;
%----- Zeeman term from transverse field:
% 1 Bohr magneton = 0.6717  Kelvin / Tesla
%                 = 0.05788 meV / Tesla
% Holmium Lande factor is 5/4.

%Hzeeman=(-5/4*0.05788*h)*Jx;

if length(h)==1
  Hzeeman=(-1.25*0.05788*h)*Jx;
else
  Hzeeman=(-1.25*0.05788)*(h(1)*Jx+h(2)*Jy+h(3)*Jz);
end

Hcfz=cf()+Hzeeman;

% Include nuclear spin through hyperfine coupling?
if 1==1
  nJ=2*8+1;
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
    hHcfz(n*nJ+(1:nJ),n*nJ+(1:nJ))=Hcfz;
    hJx(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jx;
    hJy(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jy;
    hJz(n*nJ+(1:nJ),n*nJ+(1:nJ))=Jz;
    for m=0:nI-1
      hIJ(n*nJ+(1:nJ),m*nJ+(1:nJ))=Jx*Ix(n+1,m+1)+Jy*Iy(n+1,m+1)+Jz*Iz(n+1,m+1);
    end
  end
  hIJ=0.039/11.6*hIJ;
  Hcfz=hHcfz+hIJ;
  Jx=hJx;
  Jy=hJy;
  Jz=hJz;
end;

for iterations=1:100
 % iterations 
 % Diagonalize
  h_dipol=d*[jx jy jz]';
  Ham=Hcfz-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;
  [v,e]=eig(Ham);
  e=real(diag(e));
  e=e-min(e);
  [e,n]=sort(e);
  v=v(:,n);
  % e contains eigenvalues,
  % rows of v contains eigenvectors
  if t==0
    % At zero temperature, use only lowest eigenvalue.
    jzold=jz;
    jx=real(v(:,1)'*Jx*v(:,1));
    jy=real(v(:,1)'*Jy*v(:,1));
    jz=real(v(:,1)'*Jz*v(:,1));
  else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));
    jzold=jz;
    jx=real(diag(v'*Jx*v))'*z;
    jy=real(diag(v'*Jy*v))'*z;
    jz=real(diag(v'*Jz*v))'*z;
  end
  % Itterate untill absolute change in <jz> is less than 0.0001
  if abs(jzold-jz)<0.0001
%    iterations

    return
  end
end
%iterations
return


