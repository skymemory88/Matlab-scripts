function MCscriptEQ(jobid)
%% An example.
%% Setup parameters.

infilename=sprintf('init_%d.mat',jobid);
load(infilename)
       
staggfield=0;

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
B=[[60.2258   -0.1164   -4.3280   0  -0.0019   -0.0850   -0.0227]     %Er [63   -0.55   -5.54  0.47  -0.000006   -0.1082   -0.0146]
   [-60   0.35   3.6  0  0.0004   0.07   0.006]                       %Ho
   [663   12.5   102  0  -0.62   -16   0]                             %Yb !!! Conradin's parameters !!!
   [224.3   -1.85   -11.7  -15.2  0.002   0.2645   0.1377]            %Tm
   [0 0 0 0 0 0 0]                                                    %Gd
   [0 0 0 0 0 0 0]]/1000;                                             %Y
    

F=[0.005; 0; 0; 0; 0; 0]; %cf Science LiErF4 2012, supp mat

params.staggconf=[1 0 0; -1 0 0; -1 0 0; 1 0 0];
params.staggfield=staggfield(1)*params.staggconf;
 
% Global Parameters
params.pos=[0  0  1/2 1/2      % Ions' positions within the unit cell
            0 1/2 1/2  0
            0 1/4 1/2 3/4];
params.C=[  1  1  1
           -1  1  1
           -1 -1  1
            1 -1  1];
params.replicas=[10;10;10];
params.abc=a{num_abc}; 

% Leave as it is
params.field_changed=0;

ion=cell(size(name,1), 1);

for i=1:size(name,1)
    ion{i}.name=name{i};
    ion{i}.num=i;
    ion{i}.J=J(i);
    ion{i}.L=L(i);
    ion{i}.S=S(i);
    ion{i}.B=B(i,:);  
    ion{i}.F=F(i);   
    ion{i}.gLande=gLande(ion{i}.L,ion{i}.S); 
    ion{i}.Hcf=cf(ion{i}.B,ion{i}.J);
    [ion{i}.VV,~]=Ising_basis(params.field(1,:), ion{i});
end 


%% Calculate

[lattice,params,~,lat_mom]=lat(ion,params);
disp('Lattice done.')

inter=getinter(params);
disp('Dipole interactions done.')

%% 

if(strcmp(params.meas_type,'temp'))
    disp('Simulation type: Temperature scan')
    [relaxE,EQlat_mom,lattice,params]=TloopEQ(ion,params,inter,lattice,lat_mom);
elseif(strcmp(params.meas_type,'field'))
    disp('Simulation type: Field scan')
    [relaxE,EQlat_mom,lattice,params]=FieldloopEQ(ion,params,inter,lattice,lat_mom);
end

filename=sprintf('resultsEQ_%d.mat',params.jobid);
save(filename,'inter','lattice','relaxE','params','ion','EQlat_mom');


