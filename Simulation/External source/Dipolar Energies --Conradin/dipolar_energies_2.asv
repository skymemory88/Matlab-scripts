clear all

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

rr=[rrx(:)+1/2 rry(:)+1/2 0*rrx(:)+1/2];
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
    Jj=-[cos(th(n)) sin(th(n)) 0];
%     Jj=[-1 0 0];
    E(n)=-Ji*D*Jj';
end
plot(E)
E3(5,:)=E;

h=plot(th,E3*1000)
xlabel('Angle [rad]')
ylabel('Dipole energy [arb u]')
legend(h,'FM neighbour layer','Intra-layer','AFM neighbour layer','sum',1)
%line([1.5 1.5],[-40 0],'color','k','linestyle','--')
print -dpdf dipolar_energies_layers_rotating
