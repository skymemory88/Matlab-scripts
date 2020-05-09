function create_init_file(jobid)

params.jobid=jobid;

params.field=[0 0 0];
% params.field=zeros(40,3);
% params.field(1:15,3)=linspace(0,0.8,15)';
% params.field(16:25,3)=linspace(0.84,19.8,10)';
% params.field(26:40,3)=linspace(20,24,15)';


params.temp=linspace(0,0.6,50);
% params.temp=0;
% params.temp=zeros(1,40);
% logtemp=logspace(0,-5,25);
% params.temp(1:25)=0.130*(1-logtemp);
% params.temp(26:40)=linspace(0.130,0.250,15);
% params.temp=zeros(1,40);
% params.temp(1:7)=linspace(0,0.016,7);
% params.temp(8:33)=linspace(0.020,0.040,26);
% params.temp(34:40)=linspace(0.044,0.080,7);

params.L=7;

params.prop=[1;0;0;0;0;0];
num_abc=1;
params.NiterEQ=2e6;
params.Nitermeas=1e6;

% 1 for random, 2 for ordered
params.init_Er=2;
params.init_Ho=2;
params.init_Yb=2;

% params.meas_type='field';
params.meas_type='temp';

filename=sprintf('init_%d.mat',jobid);
save(filename, 'params', 'num_abc')
end