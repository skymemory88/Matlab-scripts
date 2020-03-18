function [jx,jy,jz,chi]=remfchi(h,t,jx,jy,jz,Hcf)

% This program performs a simple meanfield calculation
% for LiHoF4


% Initiate J operators
Jz=diag(8:-1:-8);
Jp=diag(sqrt((8-[7:-1:-8]).*(8+1+[7:-1:-8])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

% Get dipole coupling J_0
%d=dipol_direct(0,10);
%d=[-0.0301383 0 0
%   0 -0.0301383 0
%   0 0 0.0602766];
%Demagnetisierung
%d=d+(4*pi/3)*0.01389*eye(3);
% Convert from AA^-3 to meV by multiplying by
% (g_Lande \mu_B)^2  (muB^2=0.05368 meV AA^3)
%d=d*(5/4)^2*0.05368;
%Lorentzfaktor 
%d=d+(4*pi/3)*0.0011654*eye(3);
% This is David Bitko's choice of coupling:
%d=2*[0 0 0;0 0 0;0 0 0.027]/11.6;

%Jensen(2000): Jex=-0.436mueV 
d=diag([3.912 3.912 6.821])/1000-4*0.436*eye(3)/1000;
%----- Zeeman term from transverse field:
% 1 Bohr magneton = 0.6717  Kelvin / Tesla
%                 = 0.05788 meV / Tesla
% Holmium Lande factor is 5/4.

xdotierung=1;
d=xdotierung*d;

%Hzeeman=(-5/4*0.05788*h)*Jx;

if length(h)==1
  Hzeeman=(-1.25*0.05788*h)*Jx;
else
  Hzeeman=(-1.25*0.05788)*(h(1)*Jx+h(2)*Jy+h(3)*Jz);
end

Hcfz=Hcf+Hzeeman;

% Include nuclear spin through hyperfine coupling?
if 1==0
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
end

for iterations=1:100
 % iterations 
  % Add mean-field
  h_dipol=d*[jx jy jz]';
  %Demagnetisierung
  %h_dipol=h_dipol-(1.25*0.05788)*4*pi*1.619*[jx jy jz]'/3/8;
  %Demagnetisierungsfaktor 8pi/3 * rel magnetisierung * 16.19 kOe
  Ham=Hcfz-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;
  % Diagonalize
  [v,e]=eig(Ham);
  e=real(diag(e));
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
    z=0*e;
    z(1)=1;
    
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

    break
  end
end
%iterations

%chi berechnen
    %k_B in meV pro kelvin%
if t==0
    beta=0;
else
    beta=11.6/t;
end

  if t>0;
      e=e-min(e);
      ne=exp(-e*11.6/t)/sum(exp(-e*11.6/t));
  else
      ne=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0  
  end

  wzaehler=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)));
  wnenner=(e*ones(1,size(v,1))-ones(size(v,1),1)*e');
  wgleich=(abs(wnenner)<0.01);
  wzaehler=wzaehler.*(wgleich~=1);
  wnenner=wnenner+wgleich;
  w=wzaehler./wnenner;
  
  
  mJx=v'*Jx*v;
  mJy=v'*Jy*v;
  mJz=v'*Jz*v; 
  
  chi0=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
        sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
        sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];
  nee=wgleich.*(ones(size(v,1),1)*ne');
  chi0=chi0+ beta*[sum(sum(mJx.'.*mJx.*nee)) sum(sum(mJx.'.*mJy.*nee)) sum(sum(mJx.'.*mJz.*nee))
                   sum(sum(mJy.'.*mJx.*nee)) sum(sum(mJy.'.*mJy.*nee)) sum(sum(mJy.'.*mJz.*nee))
                   sum(sum(mJz.'.*mJx.*nee)) sum(sum(mJz.'.*mJy.*nee)) sum(sum(mJz.'.*mJz.*nee))];
  chi0=chi0-beta*[diag(mJx)'*ne diag(mJy)'*ne diag(mJz)'*ne]'*[diag(mJx)'*ne diag(mJy)'*ne diag(mJz)'*ne];        
  
load('dipdir','dd0') %Einheit korrekt? dd0=dipole_direct([0 0 0],20)*(5/4)^2*0.05368; Lorentzfaktor integriert    d(:,:,nt,mt)=d(:,:,nt,mt)+(4*pi/3)*0.01389*eye(3)/4; 

dd0=xdotierung*dd0;

for n=1:4  %n,m Summe über Ionen in der Einheitszelle
    for m=1:4
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:)*dd0(:,:,n,m);
    end
end
chi=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:);chi0(:,:);chi0(:,:);chi0(:,:)])); 

%cxx=real(chi(1,1));
%cyy=real(chi(2,2));
%czz=real(chi(3,3));

return
if 1==1
    n=length(e);
    %k_B in meV pro kelvin%
    beta=1/11.6;
    cxx=0;
    czz=0;
    mJx=v'*Jx*v;
   % mJy=v'*Jy*v;
   mJz=v'*Jz*v;    
    for i=1:n
        for j=1:n
            if e(i)==e(j)
                cxx=cxx+beta*z(i)*mJx(i,j)*mJx(j,i); 
                czz=czz+beta*z(i)*mJz(i,j)*mJz(j,i);
            else
                cxx=cxx+mJx(i,j)*mJx(j,i)*(z(i)-z(j))/(e(j)-e(i));
                czz=czz+mJz(i,j)*mJz(j,i)*(z(i)-z(j))/(e(j)-e(i));
            end    
        end
    end
    cxx=real(cxx-beta*jx*jx);
    czz=real(czz-beta*jz*jz);
end



