%% An example.
%% Setup parameters.

global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.0; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=30; % period between exponential fit.
strategies.expfit_deltaN=5; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

tic;
Hx=[0:0.0005:0.07 0.075:0.05:2.8];
fields=[
    Hx
    zeros(size(Hx))
    zeros(size(Hx))
    ];

temp=0;

Erbium=0.0; % proportion erbium

Holmium=1.0; % proportion holmium

exEr=0; % magnetic exchange in J*S_iS_j

exHo=-0.0001; %0.1 microeV from Henrik and Jens PRB

renorm_Ho=[1,1,0.785];
renorm_Er=[1,1,1];
% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

ErHyp=0; % proportion, whithin the erbium percentage, of isotops carrying a nuclear spin (25% is default).

HoHyp=false;

demagn=true;

alpha=0; % shape is a sphere.

%% Calculate

[JEr,altJEr,JErhyp,altJErhyp,JHo,altJHo,history]=LiErxHoyY1_x_yF4(Erbium,Holmium,temp,fields,exEr,exHo,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho);

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result
figure(1)
plot(Hx,altJEr(1,:),Hx,JEr(3,:),Hx,JHo(3,:),Hx,altJHo(1,:));
legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>','<J^z_{Ho}>','<(-1)^{i}J^x_{Ho}>');
title(sprintf('LiHo_{%1.1f}Er_{%1.1f}F_4',Holmium,Erbium))
xlabel('H_x [T]');
% hgsave('benchmark1.fig')

%% Setup parameters.

fields=[0;0;0];
% 
temp=0:0.005:1.2; % temperature (or range of temperatures)

Erbium=0.6; % proportion erbium

Holmium=0.4; % proportion holmium

exEr=0; % magnetic exchange in J*S_iS_j

exHo=-0.0001; %0.1 microeV from Henrik and Jens PRB

renorm_Ho=[1,1,0.785];
renorm_Er=[1,1,1];
% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

ErHyp=0; % proportion, whithin the erbium percentage, of isotops carrying a nuclear spin (25% is default).

HoHyp=false;

demagn=true;

alpha=0; % shape is a sphere.

%% Calculate

[JEr,altJEr,JErhyp,altJErhyp,JHo,altJHo,history]=LiErxHoyY1_x_yF4(Erbium,Holmium,temp,fields,exEr,exHo,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho);

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result
figure(2)
plot(temp,altJEr(1,:),temp,JEr(3,:),temp,JHo(3,:),temp,altJHo(1,:));
legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>','<J^z_{Ho}>','<(-1)^{i}J^x_{Ho}>');
title(sprintf('LiHo_{%1.1f}Er_{%1.1f}F_4',Holmium,Erbium))
xlabel('T [K]')
% hgsave('benchmark2.fig')
fprintf('Benchmark elapsed time: %fs\n',toc)