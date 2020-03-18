%% An example.
%% Setup parameters.

global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.0; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

Hx = 0:0.1:9.0;
Hy = zeros(size(Hx));
Hz = zeros(size(Hx));
%define the components of DC magnetic fields in cartesian coordinates

fields=[Hx; Hy; Hz];
temp = 0.08; %set the temperature in kelvins
% fields=[0;0;0];
% temp=0:0.01:2.5; % set temperature range in kelvins

Erbium = 0.0; % proportion Erbium

Holmium = 1-Erbium; % proportion Holmium

exEr = 0; % magnetic exchange in J*S_iS_j

exHo = -0.0001; %0.1 ueV from Henrik and Jens PRB

renorm_Ho = [1,1,0.785]; %renormalization factor for Holmium
renorm_Er = [1,1,1]; %renormalization factor for Holmium
% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

ErHyp = 0.0; % proportion, whithin the erbium percentage, of isotops carrying a nuclear spin (25% is default).

HoHyp = true; %Hyperfine interactions for Holmium 

demagn = true; %Demagnetization factor

alpha = 0; 
% Set the sample shape for demagnetization field calculation (<1: cylinder, 1:sphere, 0: needle, inf: disk, else:'smarties')

%% Calculate

% profile on
% profile clear
[JEr,altJEr,JErhyp,altJErhyp,JHo,altJHo,history]=LiErxHoyY1_x_yF4(Erbium,Holmium,temp,fields,exEr,exHo,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho);
% profile report

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result
figure(1)
% plot(temp,altJEr(1,:),temp,JEr(3,:),temp,JHo(3,:),temp,altJHo(1,:));
plot(Hx,altJEr(1,:),Hx,JEr(3,:),Hx,JHo(3,:),Hx,altJHo(1,:));
legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>','<J^z_{Ho}>','<(-1)^{i}J^x_{Ho}>');
title(sprintf('LiHo_{%1.1f}Er_{%1.1f}F_4',Holmium,Erbium))
% xlabel('T [K]')
xlabel('H_x [T]');
print -dpng example.png

%% Plot history
cols=jet(length(temp)-1);
iterations=0;
figure(2)
for t=2:length(temp)
%    relax=squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2);
    relax=sqrt((squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2)).^2+(squeeze(history{t}.er(:,1,1))-history{t}.er(end,1,1)).^2);
    line(temp(t)*ones(size(relax)),1:size(relax,1),relax,'Color',cols(t-1,:));
    hold on
    iterations=iterations+size(relax,1);
end
fprintf('iterations=%d\n',iterations);
set(gca,'Box','on')
set(gca,'XGrid','on')
set(gca,'YGrid','on')
set(gca,'ZGrid','on')
xlabel('Temperature');
ylabel('Iterations');
zlabel('Ho J_i^z-<J^z>');
