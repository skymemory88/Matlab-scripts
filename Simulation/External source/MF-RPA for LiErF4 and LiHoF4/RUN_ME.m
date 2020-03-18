function RUN_ME
% Function to simulate the inelastic neutron scattering measurements on LiErF4 using MF-RPA
% Requires spec1d to run
% RUN_ME tested on Matlab 2012a on 16.05.2017 by PB


% Function to simulate the LET field experiment using RPA
plotopt.col  = [0.2 0.2 0.7];
plotopt.col2 = abs(0.8-plotopt.col);
plotopt.lnwd = 2;
plotopt.ftsz = 12;
plotopt.mksz = 5;


global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

% a) Define temperature
t=[0.02:0.1:1];
% b) Define energy range and DHO width
omega=[-0.1:0.001:0.25];
epsilon=0.005;

% c) Define field
% h=0:0.01:1;
h = 0.53;           % Critical field found from RPA

% d) Define Q-wavevector:
ql=[0.97];
q=ql'*[1 0 0]+ones(size(ql'))*[0 0 0];


xEr=1;
exEr=0; % magnetic exchange in J*S_iS_j
exHo=-0.0001; %0.1 microeV from Henrik and Jens PRB
renorm_Ho=[1,1,1];
renorm_Er=[1,1,1];
withdemagn=true;
alpha=1; % shape is a sphere.
%Jex=0/((6/5)^2*0.05368);
momente_ho=[0 0 1
        0 0 1
        0 0 1
         0 0 1];


for nh=1:length(t)
    hvec=[0 0 h];
    temp = t(nh);
    
    %% Domain 1 ===========================================================
    momente_er=[1 0 0
        -1 0 0
        -1 0 0
        1 0 0];
    
    momente=remf_2(hvec,temp,momente_er,momente_ho,exHo,exEr,xEr,1-xEr,0.23,0,true,alpha,renorm_Er,renorm_Ho);
    chi0=chi0_w(hvec,temp,momente,omega,epsilon,0,0.00043421);
    
        chi=chi_qw(q,chi0);
 
        for nw=1:length(omega)
            if temp==0
                Skw(nw)=1/pi* ...
                    sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
            else
                Skw(nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
                    sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
            end
        end

        srpa1(nh)=spec1d(omega,Skw,sqrt(Skw));
    
    %% Domain 2 ===========================================================
    momente_er=[0 1 0
        0 1 0
        0 -1 0
        0 -1 0];
    
    momente=remf_2(hvec,temp,momente_er,momente_ho,exHo,exEr,xEr,0,0.23,0,true,alpha,renorm_Er,renorm_Ho);
    chi0=chi0_w(hvec,temp,momente,omega,epsilon,0,0.00043421);
    
        chi=chi_qw(q,chi0);
         
        for nw=1:length(omega)
            if temp==0
                Skw(nw)=1/pi* ...
                    sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
            else
                Skw(nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
                    sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
            end
        end
        
        srpa2(nh)=spec1d(omega,Skw,sqrt(Skw));
    
    % Total contribution to the scattering cross-section from both domains
    srpa_total(nh) = srpa1(nh) + srpa2(nh);
    
end


hfig = setfig(1);

mapplot(srpa_total,t)
xlabel('$T$ (K)','fontsize',plotopt.ftsz,'interpreter','latex')
ylabel('Energy (meV)','fontsize',plotopt.ftsz,'interpreter','latex')
set(gca,'fontsize',plotopt.ftsz)
set(hfig,'colormap',mycmap)
ylim([0 0.2])


end

%% Subfunctions -----------------------------------------------------------


function hfig = setfig(nfig)

hfig = figure(nfig);
clf
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2) 600 500])

end


