%% An example.
%% Setup parameters.

clear all

global strategies; % global convergence strategies switches
strategies.powerlaw=false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

%<J> as a function of H (Q=1) or T (Q=0) ?
Q=1;
HypFr = 1.0; %Hyperfine structure scaling factor

if Q
    phi = 0*pi/180;                                 % rotation angle          
    
    Rx = [1         0           0
          0         cos(phi)    sin(phi);
          0        -sin(phi)    cos(phi)];          % rotation about x 
      
    Ry = [cos(phi)  0          -sin(phi);
          0         1           0;
          sin(phi)  0           cos(phi)];          % rotation about y 
      
    Rz = [cos(phi)  sin(phi)    0;
         -sin(phi)  cos(phi)    0;
          0         0           1];                 % rotation about z

     
%     temp=[0.08, 0.1 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.8];
    temp = 0.1;
%     Hz=[0:0.4:1.2 1.25:0.05:1.6];
    Hz=0:0.1:6;
%     Hz=2:-0.005:0;
%     Hz = [0:1:15 15.1:0.1:20 20.01:0.01:25 25:0.1:30];
    fields=Rz*[Hz
            zeros(size(Hz))
            zeros(size(Hz))
                       ];
else
    temp=0:0.1:2; %  range of temperatures
    fields=[0;0;0.1];
end



%Ions' names
ion.name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
ion.name_hyp=[{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];

%Ions' proportions
ion.prop=[0;1;0;0;0;0];

%Ions' J, L, S values
%      Er       Ho      Yb      Tm      Gd      Y
ion.J=[15/2;    8;      7/2;    6;      7/2;    1];
ion.L=[6;       6;      3;      5;      0;      1];
ion.S=[3/2;     2;      1/2;    1;      7/2;    1];

%Ions' lattice parameters
ion.a=[{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
       {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
       {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
       {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y

%Ions' cf parameters
% ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
%        [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
%        [663   12.5   102   -0.62   -16   0]                             %Yb        
%        [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
%        [0 0 0 0 0 0]                                                    %Gd
%        [0 0 0 0 0 0]]/1000;                                             %Y

% Ions' cf parameters
ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
       [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
       [663   12.5   102   -0.62   -16   0]                             %Yb
%        [646.2016   15.3409  116.4854   -0.0686  -15.1817    0.0000]     %Yb 
       [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
       [0 0 0 0 0 0]                                                    %Gd
       [0 0 0 0 0 0]]/1000;                                             %Y
 
   
%Ions' renorm
ion.renorm=[[1,1,1]         %Er
            [1,1,0.785]     %Ho
            [1,1,1]         %Yb
            [1,1,1]         %Tm
            [1,1,1]         %Gd
            [1,1,1]];       %Y

%Ions' hyperfine
ion.hyp=[0;1;0;0;0;0];

%Ions' nuclear spin I
%Ions' nuclear spin I
A_Ho=0.003361; %Henrik
% A_Ho=1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000; % 1985 paper

ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[0.00043412; A_Ho; 0; 0; 0; 0]*HypFr;

%Exchange
% ex.Er=0; % magnetic exchange in J*S_iS_j
ex.Ho=-0.0001; %0.1 microeV from Henrik and Jens PRB
% ex.Yb=0;
% ex.Tm=0;

% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?
ion.ex=[0; ex.Ho; 0; 0; 0; 0];

%Ions' initial moments
ion.mom(:,:,1)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Er
% ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1];  %Ho
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho why 2.6? (Yikai, 22.11.2019
% ion.mom(:,:,2)=[3.473 -0.045 0.1;  3.473 -0.045 0.1;  3.473 -0.045 0.1; 3.473 -0.045 0.1];  %Ho
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

demagn=true;

alpha=0; % shape is a sphere.

%% Calculate

% profile on
% profile clear
for i=1:length(temp)
    [ion,history,E]=LiIonsF4(ion,temp(i),fields,demagn,alpha);
end
% profile report

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result

    
if Q
    altJ_H=zeros([3,size(fields,2),size(ion.name,1)]);
    J_H=zeros([3,size(fields,2),size(ion.name,1)]);  
    Jnorm_H=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_H(:,:,i)=squeeze(mean(ion.altJs(:,:,:,:,i),1));
        J_H(:,:,i)=squeeze(mean(ion.Js(:,:,:,:,i),1));
        Jnorm_H(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
else
    altJ_T=zeros([3,size(temp,2),size(ion.name,1)]);
    J_T=zeros([3,size(temp,2),size(ion.name,1)]); 
    Jnorm_T=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_T(:,:,i)=mean(ion.altJs(:,:,:,i),1);
        J_T(:,:,i)=mean(ion.Js(:,:,:,i),1);
        Jnorm_T(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
end

n=find(ion.prop);
% hfig1 = figure(1);
% clf
% set(hfig1,'position',[10 10 600 600])

if Q
    p1=plot(Hz,Jnorm_H(:,1));
%     p2=plot(Hz,altJ_H(1,:,n),Hz,J_H(3,:,n));
    p2=plot(Hz,J_H(:,:,n), Hz, J_H(:,:,n));
%     legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>');
%     legend('<(-1)^{i}J^x>','<J^z>');
    legend('<J^y>','<J^z>');
    title(ion.name(n));
%     title(sprintf('LiGdF_4'))
    xlabel('H_z [T]');
    ylabel('|J|');
%     saveplots(hfig1,['Simulation of LiYbF4 Hz vs J'])
    
%     p3=plot(Hz,altJ_H(1,:,n),Hz,altJ_H(2,:,n),Hz,altJ_H(3,:,n),Hz,J_H(1,:,n),Hz,J_H(2,:,n),Hz,J_H(3,:,n));
% %     legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>');
%     legend('<(-1)^{i}J^x>','<(-1)^{i}J^y>','<(-1)^{i}J^z>','<J^x>','<J^y>','<J^z>');
%     title(char(char(strcat('Li', ion.name(n), 'F_4')),...
%             ['B = ' num2str(ion.B(3,:)*1000) ' meV'],...
%             ['H = H_0(' num2str(-sin(phi),3) ', ' num2str(cos(phi),3) ', 0), \phi(deg) = ' num2str(phi*180/pi)]))
%     
% %     title(sprintf('LiGdF_4'))
%     xlabel('H_z [T]');
else
    figure
    p=plot(temp,altJ_T(1,:,n),temp,J_T(3,:,n));
    legend('<(-1)^{i}J^x>','<J^z>');
%     title(sprintf('LiEr_{%1.1f}Ho_{%1.1f}Yb_{%1.1f}Tm_{%1.1f}Gd_{%1.1f}Y_{%1.1f}F_4',ion.prop(1),ion.prop(2),ion.prop(3),ion.prop(4),ion.prop(5),ion.prop(6)))
    title('LiYbF_4')
    xlabel('T [K]')
    ylabel('|J|');
%     saveplots(hfig1,['Simulation of LiYbF4 T vs J'])
end

% print -dpng example.png

%% Plot history
% cols=jet(length(temp)-1);
% iterations=0;
% figure(2)
% for t=2:length(temp)
% %    relax=squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2);
%     relax=sqrt((squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2)).^2+(squeeze(history{t}.er(:,1,1))-history{t}.er(end,1,1)).^2);
%     line(temp(t)*ones(size(relax)),1:size(relax,1),relax,'Color',cols(t-1,:));
%     hold on
%     iterations=iterations+size(relax,1);
% end
% fprintf('iterations=%d\n',iterations);
% set(gca,'Box','on')
% set(gca,'XGrid','on')
% set(gca,'YGrid','on')
% set(gca,'ZGrid','on')
% xlabel('Temperature');
% ylabel('Iterations');
% zlabel('Ho J_i^z-<J^z>');