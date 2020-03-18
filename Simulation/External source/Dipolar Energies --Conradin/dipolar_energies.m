% this is an analysis of how stiff spins are coupled between different
% planes of the AFM structure

clear Eml
layers=-1:1;
for nlayer=1:length(layers)
layer=layers(nlayer);

N=20;
NN=N*1.5; % slightly larger than sqrt(2)
[rrx,rry]=meshgrid(-NN:NN,-NN:NN);

Rm=1:0.5:N;
Rm=0.5:0.5:N;
Rm=0:0.1:15;

for m=1:length(Rm)
if layer==0% within same plane
Jj=[1 0 0];
rr=[rrx(:) rry(:) 0*rrx(:)];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
Emin=-19.5e-3;
end

if layer==1% AFM neigbour layer
Jj=[-1 0 0];
rr=[rrx(:) rry(:)+1/2 0*rrx(:)+1/4];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
%rr=[0 1/2 1/4
%    0 -1/2 1/4
%    1 1/2 1/4
%    -1 1/2 1/4
%    1 -1/2 1/4
%    -1 -1/2 1/4];
Emin=-23.5e-3;
end

if layer==-1 % FM neigbour layer
Jj=[1 0 0];
rr=[rrx(:)+1/2 rry(:) 0*rrx(:)-1/4];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
%rr=[1/2 0 -1/4
%    -1/2 0 -1/4
%    1/2 1 -1/4
%    -1/2 1 -1/4
%    1/2 -1 -1/4
%    -1/2 -1 -1/4];
Emin=-9.3e-3;
end

D=zeros(3,3);
for n=1:size(rr,1)
  D=D+dipole_r(rr(n,:));
end

if 0 % draw coupling versus angle
th=[0:0.01:1]*2*pi;
else
th=0;
end

E=zeros(size(th));
for n=1:length(th)
    Ji=[cos(th(n)) sin(th(n)) 0];
    E(n)=-Ji*D*Jj';
end

%clf
%plot(th,E*1000)
%Em(m)=min(E)*1000;
Eml(nlayer,m)=min(E)*1000;
end

%plot(Rm,Em,'x')

end

h=plot(Rm,Eml,'x-')
xlabel('Radius [a]')
ylabel('Dipole energy [arb u]')
legend(h,'FM neighbour layer','Intra-layer','AFM neighbour layer',1)
line([1.5 1.5],[-40 0],'color','k','linestyle','--')
print -dpdf dipolar_energies_layers

%===== try turning two layers in phase or antiphase

% result is that turning layers in pahse and antiphase between the two
% equivalent ordering patterns cost no energy ! So, the dipole coupling
% keeps full planar isotropy.
% Guess quesiton is now how about crystal field and excahnge interaction.

N=20;
NN=N*1.5; % slightly larger than sqrt(2)
[rrx,rry]=meshgrid(-NN:NN,-NN:NN);
Rm=15;m=1;

rr=[rrx(:) rry(:)+1/2 0*rrx(:)+1/4];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
D=zeros(3,3);
for n=1:size(rr,1)
  D=D+dipole_r(rr(n,:));
end
th=[0:0.01:1]*2*pi;
E=zeros(size(th));
for n=1:length(th)
    Ji=[cos(th(n)) sin(th(n)) 0];
    Jj=[-cos(th(n)) sin(th(n)) 0];
%     Jj=[-1 0 0];
    E(n)=-Ji*D*Jj';
end
E3(1,:)=E;

rr=[rrx(:) rry(:) 0*rrx(:)];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
D=zeros(3,3);
for n=1:size(rr,1)
  D=D+dipole_r(rr(n,:));
end
th=[0:0.01:1]*2*pi;
E=zeros(size(th));
for n=1:length(th)
    Ji=[cos(th(n)) sin(th(n)) 0];
    Jj=[cos(th(n)) sin(th(n)) 0];
%    Jj=[1 0 0];
    E(n)=-Ji*D*Jj';
end
E3(2,:)=E;

rr=[rrx(:)+1/2 rry(:) 0*rrx(:)-1/4];
rr(sum(rr.^2,2)==0,:)=[];
rr(sum(rr.^2,2)>Rm(m)^2,:)=[];
D=zeros(3,3);
for n=1:size(rr,1)
  D=D+dipole_r(rr(n,:));
end
th=[0:0.01:1]*2*pi;
E=zeros(size(th));
for n=1:length(th)
    Ji=[cos(th(n)) sin(th(n)) 0];
    Jj=[cos(th(n)) -sin(th(n)) 0];
%    Jj=[1 0 0];
    E(n)=-Ji*D*Jj';
end
E3(3,:)=E;

E3(4,:)=sum(E3,1);

plot(th,E3*1000)
xlabel('Angle [rad]')
ylabel('Dipole energy [arb u]')
legend(h,'FM neighbour layer','Intra-layer','AFM neighbour layer','sum',1)
%line([1.5 1.5],[-40 0],'color','k','linestyle','--')
print -dpdf dipolar_energies_layers_rotating
