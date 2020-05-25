clearvars
alpha=linspace(0,pi,201);
beta=linspace(0,2*pi,60);

mom=zeros(3,length(alpha), length(beta));
E=zeros(length(alpha), length(beta));

field=[0,0,0];

%Ions' names
name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];

%Ions' J, L, S values
J=[15/2; 8; 7/2; 6; 7/2; 1];
L=[6; 6; 3; 5; 0; 1];
S=[3/2; 2; 1/2; 1; 7/2; 1];

%Ions' lattice parameters
a=[{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
   {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
   {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
   {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
   {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
   {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y

%Ions' cf parameters
B=[[60.2258   -0.1164   -4.3280   0  -0.0019   -0.0850   -0.0227]     %Er [60.2258   -0.1164   -4.3280   0  -0.0019   -0.0850   -0.0227]  [63   -0.55   -5.54  0.47  -0.000006   -0.1082   -0.0146]
   [-60   0.35   3.6  0  0.0004   0.07   0.006]                       %Ho
   [663   12.5   102  0  -0.62   -16   0]                             %Yb !!! Conradin's parameters !!!
   [224.3   -1.85   -11.7  -15.2  0.002   0.2645   0.1377]            %Tm
   [0 0 0 0 0 0 0]                                                    %Gd
   [0 0 0 0 0 0 0]]/1000;                                             %Y


for i=2 % choose the ion
    ion{i}.name=name{i};
    ion{i}.num=i;
    ion{i}.J=J(i);
    ion{i}.L=L(i);
    ion{i}.S=S(i);
    ion{i}.B=B(i,:);
    ion{i}.gLande=gLande(ion{i}.L,ion{i}.S);
    ion{i}.Hcf=cf(B(i,:),J(i));
    [ion{i}.VV,~]=Ising_basis([0,0,0],ion{i});
end 

params.field_changed=1;

for j=1:length(alpha)
    for k=1:length(beta)
    [mom(:,j,k),E(j,k),params]=cl_eff_mod(2*alpha(j), beta(k), field, ion{i}, params);
    end
end

Y=zeros(length(alpha)*length(beta),1);
X=zeros(length(alpha)*length(beta),1);
Z=zeros(length(alpha)*length(beta),3);
theta2=zeros(length(alpha)*length(beta),1);
Ecf=zeros(length(alpha)*length(beta),1);
theta=zeros(length(alpha),length(beta));
phi=zeros(length(alpha),length(beta));
normes=zeros(length(alpha),length(beta));

for j=1:length(alpha)
    for k=1:length(beta)
         X((j-1)*60+k)=alpha(j);
         Y((j-1)*60+k)=beta(k);
           
         Z((j-1)*60+k,1)=mom(1,j,k);
         Z((j-1)*60+k,2)=mom(2,j,k);
         Z((j-1)*60+k,3)=mom(3,j,k);
        [phi(j,k), theta(j,k),normes(j,k)]=cart2sph(mom(1,j,k),mom(2,j,k),mom(3,j,k));
         theta2((j-1)*60+k)=theta(j,k);
         Ecf((j-1)*60+k)=E(j,k);
    end
end
%%
b_ind=1;

figure(1)
plot(alpha,squeeze(mom(1,:,b_ind)),'-s','Markersize',2)
xlabel('\alpha')
ylabel('Jx')

figure(2)
plot(alpha,squeeze(mom(2,:,b_ind)),'-s','Markersize',2)
xlabel('\alpha')
ylabel('Jy')

figure(3)
plot(alpha,squeeze(mom(3,:,b_ind)),'-s','Markersize',2)
xlabel('\alpha')
ylabel('Jz')

% b_ind=1;

figure(4)
plot(alpha,theta(:,b_ind),'-s','Markersize',2)
xlabel('\alpha')
ylabel('\theta')
tit=sprintf('beta = %f',beta(b_ind));
axis([0,pi,-pi/2,pi/2])
grid on
title(tit)

figure(5)
plot(beta,phi(b_ind,:),'-s','Markersize',2)
xlabel('\beta')
ylabel('\phi')
tit=sprintf('alpha = %f',alpha(b_ind));
axis([0,2*pi,-pi,pi])
grid on
title(tit)

figure(6)
plot3(X,Y,squeeze(Z(:,1)), 'x', 'Markersize',2)
xlabel('\alpha')
ylabel('\beta')
zlabel('J_x')
title(['Li', ion{1,i}.name,'F_4'])
grid on

figure(7)
plot3(X,Y,squeeze(Z(:,2)), 'x', 'Markersize',2)
xlabel('\alpha')
ylabel('\beta')
zlabel('J_y')
title(['Li', ion{1,i}.name,'F_4'])
grid on


figure(8)
plot3(X,Y,squeeze(Z(:,3)), 'x', 'Markersize',2)
xlabel('\alpha')
ylabel('\beta')
zlabel('J_z')
title(['Li', ion{1,i}.name,'F_4'])
grid on

figure(9)
plot3(X,Y,Ecf, 'x', 'Markersize',2)
xlabel('\alpha')
ylabel('\beta')
zlabel('E_{cf}')
title(['Li', ion{1,i}.name,'F_4'])
grid on

 figure(10)
 plot3(squeeze(Z(:,1)),squeeze(Z(:,2)),squeeze(Z(:,3)),'x', 'Markersize',2)
 title(['Li', ion{1,i}.name,'F_4'])
 grid on

%
 figure(11)
 plot3(X,Y,theta2,'x','Markersize',2)
 xlabel('\alpha')
 ylabel('\beta')
 zlabel('\theta')

% colnames='alpha theta';
% fid=fopen('theta_cl_eff_mod_Yb.txt','wt');
% fprintf(fid, '%s \n',colnames);
% for j=1:201
%     fprintf(fid, '%f %f \n', alpha(j), theta(j,b_ind));
% end
% fclose(fid);

% colnames='Jx Jy Jz';
% fid=fopen('mom3d_cl_eff_mod_Ho.txt','wt');
% fprintf(fid, '%s \n',colnames);
% for j=1:length(Z)
%     fprintf(fid,'%e %e %e \n', Z(j,1),Z(j,2),Z(j,3));
% end
% fclose(fid);
