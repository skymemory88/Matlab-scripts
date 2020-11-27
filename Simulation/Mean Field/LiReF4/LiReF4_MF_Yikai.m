function LiReF4_MF_Yikai
clearvars
%% Setup parameters.
global strategies; % global convergence strategies switches
strategies.powerlaw=false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.05; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=30; % period between exponential fit.
strategies.expfit_deltaN=5; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.
strategies.symmetry = true; % Copy the state of one site to its crystal symmetry equivalent sites instead of calculating them individually

global rundipole % Run dipole calculations or not
rundipole = true;
%% define temp / field scan
global Options;
Options.scantype = true; %Choose scan type: field scan (TRUE) or temperature scan (FALSE)
Options.hyperfine = true; % Including/excluding hyperfine calculations
Options.plotting = false; %Choose whether or not to plot the figures at the end
Options.saving = true ;

if Options.scantype == true % Temperature scan
%         temp = [0.1 0.12 0.15 0.2 0.24 0.3 0.35 0.4 0.45 0.5 0.8 1.2 1.6 1.8 2];
        temp =  0.1;
        hypFr = 1.0; % Scaling factor for hyperfine interaction
        Hmin = 0.0; % Minimum magnetic field
        Hmax = 16.0; % Maximum magnetic field
        H_step = 0.02; % Field scan resolution
        theta = 0.0/180*pi; % theta(in radian) = 0 indicates a transverse magnetic field
        phi = 60.0/180*pi; % ab-plane rotation, phi(in radian) = 0 means H along x
        
        % Calcualte each component of the DC magnetic field
        Hx = (Hmin:H_step:Hmax)*cos(phi)*cos(theta);
        Hy = (Hmin:H_step:Hmax)*sin(phi)*cos(theta);
        Hz = (Hmin:H_step:Hmax)*sin(theta);
        fields = [Hx;Hy;Hz];
elseif Options.scantype == false % Field scan
        hypFr = 1.0; % Scaling factor for hyperfine interaction
        H = 0; % Static external magnetic field
        theta = 0.0/180*pi; % theta(in radian) = 0 indicates a transverse magnetic field
        phi = 0.0/180*pi; % ab-plane rotation, phi(in radian) = 0 means H along x
        Hx = H*cos(phi)*cos(theta);
        Hy = H*sin(phi)*cos(theta);
        Hz = H*sin(theta);
        temp = 0:0.04:2.0 ; %  range of temperatures
        fields = [Hx;Hy;Hz];
else
    error('Options.scantype must equal to 0 or 1')
end
%% define LiRF4

%Ions' names
ion.name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
ion.name_hyp=[{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];

%Ions' proportions
ion.prop=[0;1;0;0;0;0]; % LiHoF4
% ion.prop=[1;0;0;0;0;0]; % LiErF4

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

% Ions' cf parameters
ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
       [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
       [646.2016   15.3409  116.4854   -0.0686  -15.1817    0.0000]     %Yb    
%        [663   12.5   102   -0.62   -16   0]                             %Yb
       [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
       [0 0 0 0 0 0]                                                    %Gd
       [0 0 0 0 0 0]]/1000.0;                                           %Y
 
%Ions' renorm
ion.renorm=[[1,1,1]         %Er
            [1,1,0.785]     %Ho
            [1,1,1]         %Yb
            [1,1,1]         %Tm
            [1,1,1]         %Gd
            [1,1,1]];       %Y


%Ions' nuclear spin I
A_Ho=0.003361; %Henrik
% A_Ho=1.6*1e9*6.62606957e-34/(1.602176565e-19)*1000; % 1985 paper
A_Er=0.00043412;

ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[A_Er; A_Ho; 0; 0; 0; 0]*hypFr;

%Exchange
ex.Er=0; % magnetic exchange in J*S_iS_j
ex.Ho=-0.0001; %0.1 microeV from Henrik and Jens PRB
% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?
ex.Yb=0;
ex.Tm=0;

% ion.ex=[0; ex.Ho; 0; 0; 0; 0]; % LiHoF4
% ion.ex=[ex.Er; 0; 0; 0; 0; 0]; % LiErF4
ion.ex=[ex.Er; ex.Ho; ex.Yb; ex.Tm; 0; 0];

%Ions' initial moments
ion.mom(:,:,1)=[1 0 0; 1 0 0; -1 0 0; -1 0 0];      %Er
% ion.mom(:,:,1)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Er from QMC code (MCscriptEQ.m)
% ion.mom(:,:,1)=[0 1 0; 0 1 0; 0 -1 0; 0 -1 0];      %Er
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho %Why multiply 2.6 (Yikai)?
% ion.mom(:,:,2)=[3.473 -0.045 0.1;  3.473 -0.045 0.1;  3.473 -0.045 0.1; 3.473 -0.045 0.1];  %Ho
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

%Ions' hyperfine
if Options.hyperfine
    % ion.hyp=[0;1;0;0;0;0]; % LiHoF4
    % ion.hyp=[1;0;0;0;0;0]; % LiErF4
    ion.hyp=[1;1;0;0;0;0]; % Er, Ho, Yb, Tm, Gd, Y
else
    ion.hyp = [0;0;0;0;0;0];
end

demagn=true; % Demagnetization factor included if TRUE
alpha=0; % shape in calculating demagnetization factor is a needle (0), sphere (1)
%% Calculate and save results inside the LiIonsF4 function
[ion,~,E,V]=LiIonsF4(ion,temp,fields,phi,theta,demagn,alpha);
%% Plot data
n = ion.prop~=0;
if Options.plotting == true
    if Options.scantype == true % For field scan
        J_H=zeros([3,size(fields,2),size(ion.name,1)]);
        Jnorm_H=zeros([size(fields,2),size(ion.name,1)]);
    else % for temperature scan
        J_H=zeros([3,size(temp,2),size(ion.name,1)]);
        Jnorm_H=zeros([size(temp,2),size(ion.name,1)]);
    end
    
    for ii=1:size(ion.name,1)
        J_H(:,:,ii)=squeeze(mean(ion.Js_hyp(:,:,:,ii),1)); % With hyperfine
        %         J_H(:,:,i)=squeeze(mean(ion.Js(:,:,:,:,i),1)); % Without hyperfine
        Jnorm_H(:,ii)=squeeze(mean(ion.Jmom_hyp_norm(:,ii),1)); % with hyperfine
        %         Jnorm_H(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1)); % Without hyperfine
    end
    
    HH = zeros(1,size(fields,2)); % Calculate the norm of the applied magnetic field
    for ii = 1:size(fields,2)
        HH(ii) = norm(fields(:,ii));
    end
    
    %Plot spin moments
    if Options.scantype == true % For field scan 
        figure
        hold on
        p1 = plot(HH(1:10:end),Jnorm_H(1:10:end,n),'-o','Color', 'black', 'LineWidth', 2, 'MarkerSize',4);
        title(sprintf('Norm of spin moments for %s', char(ion.name(n))));
        xlabel('Magnetic field (|B| [T])');
        ylabel('Spin moment norm <J_n>');
        
        figure
        hold on
        p2_1 = plot(HH,J_H(:,:,n),'-','LineWidth', 2);
        p2_2 = plot(HH(1:10:end), J_H(:,1:10:end,n),'o','MarkerSize',4); % Add spaced markers to the <J> curves
        p2_3 = plot(HH(1:10:end), Jnorm_H(1:10:end,n),'-s','MarkerSize',4);
        title(sprintf('Expectation values of spin moments for %s', char(ion.name(n))));
        legend('<J_x>','<J_y>','<J_z>','<J_{norm}>');
        xlabel('Magnetic field (|B| [T])');
        ylabel('Spin moment <J>');
    else
        figure
        hold on
        p1 = plot(temp,ion.altJmom(:,:,n),'-o','Color', 'black', 'LineWidth', 2, 'MarkerSize',4);
        title(sprintf('Alternating spin moments for %s', char(ion.name(n))));
        xlabel('Temperature (K)');
        ylabel('Alternating Spin moment <altJ_n>');
        
        figure
        hold on
        p2_1 = plot(temp,J_H(:,:,n),'-','LineWidth', 2);
        p2_2 = plot(temp,J_H(:,:,n),'o','MarkerSize',4); % Add spaced markers to the <J> curves
        p2_3 = plot(temp(1:10:end),Jnorm_H(1:10:end,n),'-s','MarkerSize',4);
        title(sprintf('Expectation values of spin moments for %s', char(ion.name(n))));
        legend('<J_x>','<J_y>','<J_z>','<J_{norm}>');
        xlabel('Magnetic field (|B| [T])');
        ylabel('Spin moment <J>');
    end
end
%% Save the data files
% if more than one temperatures provided, the data is saved insied LiIonF4() function
if Options.saving == true
    eee = squeeze(E);
    fff = fields;
    ttt = temp;
    vvv = squeeze(V);
    if Options.scantype == 1 && length(temp) == 1 % Field scan: Options.scantype=1, Temperature scan: Options.scantype=0
        elem = [ion.name(n)];
        tit = strcat('Hscan_','Li',elem(1:end),sprintf('F4_%1$3.3fK_%2$.1fDeg_%3$.1fDeg.mat', temp, theta*180/pi, phi*180/pi));
%         cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output')
        filepath = 'My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I';
        % if more than one field value provided, the data is saved insied LiIonF4() function
    elseif Options.scantype == 0 && size(fields,2) == 1 % Field scan: Options.scantype=1, Temperature scan: Options.scantype=0
        filepath = 'My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I';
        elem = [ion.name(n)];
        tit = strcat('Tscan_Li',elem(1:end),sprintf('F4_%1$3.3fK_%2$.1fDeg_%3$.1fDeg.mat', temp, theta*180/pi, phi*180/pi));
    end
    save(fullfile(filepath,char(tit)),'ttt','fff','eee','vvv','ion')
    cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4')
end
end