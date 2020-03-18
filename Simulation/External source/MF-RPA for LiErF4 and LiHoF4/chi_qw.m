function chi=chi_qw(q,chi0,Jex)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.  

if nargin<3
    Jex=0;
end

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,size(q,1));

% load('dipolwerte', 'Dpq', 'qvektor') % load file dipolwerte.mat, which contain variables Dpq and qvektor

%q=qtemp;
for n=1:size(q,1)
     D(:,:,:,:,n)=dipole_direct(q(n,:),50)+exchange(q(n,:),Jex);
end

%clear 'Dpq' 'qvektor'

% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*(6/5)^2*0.05368;%*vol;

chi=zeros(3,3,size(chi0,4),size(q,1));
for nw=1:size(chi0,4) 
for nq=1:size(q,1)
  M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
  for n=1:N  %n,m Summe �ber Ionen in der Einheitszelle
    for m=1:N
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,n,nw)*D(:,:,n,m,nq);
    end
  end
  chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:,1,nw);chi0(:,:,2,nw);chi0(:,:,3,nw);chi0(:,:,4,nw)]));
end
end


return
