function chimacroscopisch()
h=0;
t=1.8;
epsilon=0.1;
Hcf=cf();
omega=[0.01:0.01:1]; %1 meV = 0.2418 Thz
%temp=[0.1:0.1:5]
temp=t*ones(size(omega));
%omega=100*ones(size(temp))*6.5820e-013;
[e,v]=mfieldz([h 0 0],t,1,0,1,Hcf,0);
chi0=zeros(3,3,length(omega));

for n=1:length(omega)
  chi0(:,:,n)=CHI0_W(h,temp(n),omega(n),epsilon);
end

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N);

D(:,:,:,:)=dipole([0 0 0],15);

% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
vol=287.8917;
D=D*(5/4)^2*0.05368;%*vol;

chi=zeros(3,3,length(omega));
for nw=1:length(omega) 
  M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
  for n=1:N  %n,m Summe über Ionen in der Einheitszelle
    for m=1:N
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,nw)*D(:,:,n,m);
    end
  end
  chi(:,:,nw)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw)]));
end
%omega=omega/(6.5820e-013);
plot(omega,squeeze(real(chi(3,3,:))),omega,squeeze(imag(chi(3,3,:))))
figure
plot(omega,squeeze(real(chi0(3,3,:))),omega,squeeze(imag(chi0(3,3,:))))
