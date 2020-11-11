function create_init_file(jobid, LatticeSize)

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

params.temp(1:4)=linspace(0.05,0.1,4);
params.temp(5:28)=linspace(0.11,0.15,24);
params.temp(29:32)=linspace(0.16,0.25,4);
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
params.L = LatticeSize;

params.prop=[1;0;0;0;0;0]; % {'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}
num_abc=1; % {'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}

params.NiterEQ = 5e4; % Thermalization steps (in unit of lattice size)
% params.Nitermeas = 2e2; % Sampling steps (in unit of lattice size)
% params.pt_intv = 100; % Interval between parallel temperature trials
% params.meas_intv = 1000; % Interval between measurements

% for debugging
% params.NiterEQ = 5; % Thermalization steps (scaled according to the system size)
% params.Nitermeas = 5; % Sampling stemps
% params.pt_intv = 1; % Interval between parallel temperature trials
% params.meas_intv = 1; % Interval between measurements

% 1 for random, 2 for ordered, called in lat()
params.init_Er=2;
% params.init_Ho=2;
% params.init_Yb=2;

filename=sprintf(['init_',jobid,'.mat']);
save(filename, 'params', 'num_abc')
end