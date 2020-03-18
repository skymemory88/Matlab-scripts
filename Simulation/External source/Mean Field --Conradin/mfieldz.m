% identisch mit remf, ruckgabewerte sind MF-Zustaende (EW und EV) 
% zur RPA berechnung von chi_0W. Wichtig Version von remf und mfieldz muessen
% uebereinstimmen

function [e,v]=mfieldz(h,t,jx,jy,jz,Hcf,bekannt)

% This program performs a simple meanfield calculation
% for LiHoF4

%pruefen ob h,t bereits in der Liste 

if 1==0
if nargin>6
    if length(h)==1
        h=[h 0 0];
    end
    load('werte_ohnehyp', 'temperatur', 'punkte', 'Swerte')
    ind = find(t==temperatur(:));
    if ~isempty(ind)
        ind2=find((punkte(ind(1),:,1)==h(1))&(punkte(ind(1),:,2)==h(2))&(punkte(ind(1),:,3)==h(3)));
        if ~isempty(ind2)
            jx=Swerte(ind(1),ind2(1),1);
            jy=Swerte(ind(1),ind2(1),2);
            jz=Swerte(ind(1),ind2(1),3);
            bekannt=1;
        else
            bekannt=0;
        end
    else
        bekannt=0;
    end
else
    bekannt=0;
end
else
    bekannt=0;
end

% Initiate J operators
Jz=diag(8:-1:-8);
Jp=diag(sqrt((8-[7:-1:-8]).*(8+1+[7:-1:-8])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

% Get dipole coupling J_0
%d=dipol_direct(0,10);
%d=[-0.0301383 0 0
 %  0 -0.0301383 0
  % 0 0 0.0602766];

% Convert from AA^-3 to meV by multiplying by
% (g_Lande \mu_B)^2  (muB^2=0.05368 meV AA^3)
%d=d*(5/4)^2*0.05368;
% This is David Bitko's choice of coupling:
%d=2*[0 0 0;0 0 0;0 0 0.027]/11.6;
%Jensen(2000):
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

if bekannt==0 % meanfield muss noch berechnet werden
    
for iterations=1:100
 % iterations 
  % Add mean-field
  h_dipol=d*[jx jy jz]';
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

else %Meanfield aus Liste
  bekannt
  h_dipol=d*[jx jy jz]';
  Ham=Hcfz-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;
  % Diagonalize
  [v,e]=eig(Ham);
  e=real(diag(e));
  e=e-min(e);
  [e,n]=sort(e);
  v=v(:,n);    

end    
return
