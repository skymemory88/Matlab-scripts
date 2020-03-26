function ik_LiHoF4_energy_absorption
clearvars

%% Setup parameters.
global strategies; % global convergence strategies switches
strategies.powerlaw=false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.01; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=30; % period between exponential fit.
strategies.expfit_deltaN=5; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

global rundipole
rundipole = true;

%% define temp / field scan

Q=1; %<J> as a function of H (Q=1) or T (Q=0) ?

if Q == 1
        temp = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.8 1.2];
%         temp = 0.300;
        hypFr = 1.0; % Scaling factor for hyperfine interaction
        Hmin = 0.0; % Minimum magnetic field
        Hmax = 17.0; % Maximum magnetic field
        phi = 0; % phi=0 means H along x (in radian)
        theta = 0; % theta indicates a transverse magnetic field

        % Calcualte each component of the DC magnetic field
        Hx = (Hmin:0.05:Hmax)*cos(phi)*cos(theta);
        Hy = (Hmin:0.05:Hmax)*sin(phi)*cos(theta);
        Hz = (Hmin:0.05:Hmax)*sin(theta);
        
        fields = [Hx; Hy; Hz];

else
        Hx = 4;
        Hy = 0;
        Hz = 0;
        temp=0:0.1:4; %  range of temperatures
        fields=[Hx;Hy;Hz];
end


%% define LiRF4

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
       [646.2016   15.3409  116.4854   -0.0686  -15.1817    0.0000]     %Yb        
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
A_Ho=0.003361; %Henrik
% A_Ho=1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000; % 1985 paper

ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[0.00043412; A_Ho; 0; 0; 0; 0]*hypFr;

%Exchange
% ex.Er=0; % magnetic exchange in J*S_iS_j
ex.Ho = -0.0001; %0.1 microeV from Henrik and Jens PRB
% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?
% ex.Yb=0;
% ex.Tm=0;

ion.ex=[0; ex.Ho; 0; 0; 0; 0];

%Ions' initial moments
ion.mom(:,:,1)=[1 0 0; 1 0 0; -1 0 0; -1 0 0];      %Er
% ion.mom(:,:,1)=[0 1 0; 0 1 0; 0 -1 0; 0 -1 0];      %Er
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho %Why multiply 2.6 (Yikai)?
% ion.mom(:,:,2)=[3.473 -0.045 0.1;  3.473 -0.045 0.1;  3.473 -0.045 0.1; 3.473 -0.045 0.1];  %Ho
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

ion.mom_hyp=ion.mom;

demagn=true;
% alpha=1; % shape is a sphere.
alpha=0; % shape is a needle. IT WONT CONVERGE FOR SPHERE!!!


%% Calculate and save results inside the LiIonsF4 function
[ion,history,E,V]=LiIonsF4(ion,temp,fields,demagn,alpha);
%% Plot data
n = find(ion.prop~=0);
J_H=zeros([3,size(fields,2),size(ion.name,1)]);  
Jnorm_H=zeros([size(fields,2),size(ion.name,1)]);
for i=1:size(ion.name,1)
    J_H(:,:,i)=squeeze(mean(ion.Js_hyp(:,:,:,:,i),1)); % With hyperfine
%         J_H(:,:,i)=squeeze(mean(ion.Js(:,:,:,:,i),1)); % Without
%         hyperfine 

    Jnorm_H(:,i)=squeeze(mean(ion.Jmom_hyp_norm(:,:,i),1)); % with hyperfine
%         Jnorm_H(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1)); % Without hyperfine
end

HH = zeros(1,size(fields,2)); % Calculate the norm of the applied magnetic field
for ii = 1:size(fields,2)
    HH(ii) = norm(fields(:,ii));
end

% %% %Plot spin moments
% figure
% hold on
% p1 = plot(HH(1:10:end),Jnorm_H(1:10:end,n),'-o','Color', 'black', 'LineWidth', 2, 'MarkerSize',4);
% title(sprintf('Norm of spin moments for %s', char(ion.name(n))));
% xlabel('Magnetic field (|B| [T])');
% ylabel('Spin moment norm <J_n>');
% 
% figure
% hold on
% p2_1 = plot(HH,J_H(:,:,n),'-','LineWidth', 4);
% p2_2 = plot(HH(1:10:end), J_H(:,1:10:end,n),'o','MarkerSize',8); % Add spaced markers to the <J> curves
% title(sprintf('Expectation values of spin moments for %s', char(ion.name(n))));
% legend('<J_{norm}>','<J_x>','<J_y>','<J_z>');
% xlabel('Magnetic field (|B| [T])');
% ylabel('Spin moment <J>');

%%
% Save additional results
if Q == 1
    elem = [ion.name(n)];  
    tit = strcat('fieldscan_',num2str(hypFr*100,3),'%_Li',elem(1:end),'F4.mat');
    cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output')
    save(char(tit),'fields','ion','history','E','V')
else
    cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output')
    elem = [ion.name(n)];  
    tit = strcat('tempscan_',num2str(hypFr*100,3),'%_Li',elem(1:end),'F4.mat');
    save(char(tit),'temp','fields','ion','history','E','V')
end
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan')
end