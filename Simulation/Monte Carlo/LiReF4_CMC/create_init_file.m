function create_init_file(jobid)

params.jobid=jobid;

% Select simulation type and set parameters
% params.meas_type='field';
params.meas_type='temp';

params.field=[0 0 0];
% params.field=zeros(30,3);
% params.field(:,3)=linspace(0,0.9,30);
% params.field(1:15,3)=linspace(0,0.8,15)';
% params.field(16:25,3)=linspace(0.84,19.8,10)';
% params.field(26:40,3)=linspace(20,24,15)';

% params.temp(1:6) = linspace(0.1,0.6,6); % for debugging

params.temp(1:5)=linspace(0.05,0.1,5);
params.temp(6:27)=linspace(0.11,0.15,22);
params.temp(28:32)=linspace(0.16,0.25,5);
% params.temp=0.01;
% params.temp=zeros(1,40);
% logtemp=logspace(0,-5,25);
% params.temp(1:25)=0.130*(1-logtemp);
% params.temp(26:40)=linspace(0.130,0.250,15);
% params.temp=zeros(1,40);
% params.temp(1:7)=linspace(0,0.016,7);
% params.temp(8:33)=linspace(0.020,0.040,26);
% params.temp(34:40)=linspace(0.044,0.080,7);

% Supercell size
params.L=15;

params.prop=[1;0;0;0;0;0]; % {'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}
num_abc=1; % {'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}

params.NiterEQ = 5e2*4*params.L^3; % Thermalization steps (scaled according to the system size)
params.Nitermeas = 2e3*4*params.L^3; % Sampling stemps

% % for debugging
% params.NiterEQ = 4*params.L^3; % Thermalization steps (scaled according to the system size)
% params.Nitermeas = 10*4*params.L^3; % Sampling stemps

% 1 for random, 2 for ordered, called in lat()
params.init_Er=2;
% params.init_Ho=2;
% params.init_Yb=2;

filename=sprintf(['init_',jobid,'.mat']);
save(filename, 'params', 'num_abc')
end