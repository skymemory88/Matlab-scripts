function ik_LiHoF4_mf_sim_2D
clearvars -global
cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature')

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

Q = 2; %<J> as a function of H (Q=1) or T (Q=0) ?

if Q == 0  % tempscan
    
        temp=0:0.1:5; %  range of temperatures
        fields=[3.5;0;0];
        
        
else if Q == 1 % fieldscan
        
        temp=0.15;

        Hmin= 0;
        Hmax= 9;
        phi=0; % phi=0 is H along x 
        phi = phi*pi/180;
        
        Hx = [Hmin:0.1:Hmax]*cos(phi);
        Hz = [Hmin:0.1:Hmax]*sin(phi);
        Hy = 0.*Hx;
        
        fields = [Hx; Hy; Hz];

        
else if Q == 2 % 2D field-temp-scan

        temp=0.00:0.05:2.5;

        Hmin= 0;
        Hmax= 9;
        phi=0; % phi=0 is H along x 
        phi = phi*pi/180;
        
        Hx = [Hmin:0.01:Hmax]*cos(phi);
        Hz = [Hmin:0.01:Hmax]*sin(phi);
        Hy = 0.*Hx;
        
        fields = [Hx; Hy; Hz];
        
    end     
    end
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
ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[0.00043412; 0.003361; 0; 0; 0; 0];

% A_Ho_Henrik=0.003361;
% A_Ho_1985=1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000;
% A_Ho_1985=0.0066;

%Exchange
ion.ex=[0; -0.0001; 0; 0; 0; 0];
% ex.Er=0; % magnetic exchange in J*S_iS_j
% ex.Ho=-0.0001; %0.1 microeV from Henrik and Jens PRB
% ex.Yb=0;
% ex.Tm=0;

% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

%Ions' initial moments
ion.mom(:,:,1)=[1 0 0; 1 0 0; -1 0 0; -1 0 0];      %Er
% ion.mom(:,:,1)=[0 1 0; 0 1 0; 0 -1 0; 0 -1 0];      %Er
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho
% ion.mom(:,:,2)=[3.473 -0.045 0.1;  3.473 -0.045 0.1;  3.473 -0.045 0.1; 3.473 -0.045 0.1];  %Ho
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

ion.mom_hyp=ion.mom;

demagn=true;

% alpha=1; % shape is a sphere.
alpha=0; % shape is a needle. IT WONT CONVERGE FOR SPHERE!!!


%% Calculate
[ion,history,E,V]=LiIonsF4(ion,temp,fields,demagn,alpha);


%% Save results
if Q == 1
savefile = ['C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\fieldscan_',...
    num2str(temp),'_K.mat'];
    save(savefile,'temp','fields','ion','history','E','V')
elseif Q == 0
savefile = ['C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\tempscan_',...
    'Hx=',num2str(fields(1,1)),',Hy=',num2str(fields(2,1)),',Hz=',num2str(fields(3,1)),'_T.mat'];
    save(savefile,'temp','fields','ion','history','E','V')
elseif Q == 2
savefile = ['C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\field-temp-scan-10092014.mat'];
%     save(savefile,'temp','fields','ion','history','E')
    save(savefile,'temp','fields','ion','E','V')
end
cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature')


end