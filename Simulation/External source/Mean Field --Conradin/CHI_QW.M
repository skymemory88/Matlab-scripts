function chi=chi_qw(q,chi0)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.  

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,size(q,1));

load('dipolwerte', 'Dpq', 'qvektor')

%q=qtemp;
for n=1:size(q,1)
  mxyz=find((q(n,1)==qvektor(:,1)) & (q(n,2)==qvektor(:,2)) &  (q(n,3)==qvektor(:,3)));
   if ~isempty(mxyz)
     D(:,:,:,:,n)=Dpq(:,:,:,:,mxyz(1));
   else
     q(n,:)
     warning('q nicht in der Liste')
     D(:,:,:,:,n)=dipole(q(n,:),12);
     dipolq(q(n,:));
  end
end

%clear 'Dpq' 'qvektor'

% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*(5/4)^2*0.05368;%*vol;

chi=zeros(3,3,size(chi0,3),size(q,1));
for nw=1:size(chi0,3) 
for nq=1:size(q,1)
  M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
  for n=1:N  %n,m Summe über Ionen in der Einheitszelle
    for m=1:N
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,nw)*D(:,:,n,m,nq);
    end
  end
  chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw)]));
end
end


return
