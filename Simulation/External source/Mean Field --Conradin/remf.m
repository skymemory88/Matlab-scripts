function [jx,jy,jz]=remf(h,t,jx,jy,jz)

% This program performs a simple meanfield calculation
% for LiHoF4
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


for iterations=1:100
 % iterations 
 % Diagonalize
  Ham=H3Level([h,0,0],[jx,jy,jz]);
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


