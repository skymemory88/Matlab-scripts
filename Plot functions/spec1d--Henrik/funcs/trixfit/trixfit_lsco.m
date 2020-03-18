%
% TRIXFIT script
%
% Fitting triple-axis data, including resolution effects calculated using the
% Cooper-Nathans or Popovici methods.
%
% Notes: 
% (1)  The Monte Carlo integration routine used for the convolution introduces
%      noise, and as a result it is necessary to balance the parameters used
%      by the fit routine, the three fcp parameters:
%                [parameter_step no_of_iterations convergence_criterion]
%      against the number of sampling points used by the Monte Carlo. In general,
%      the larger the number of sampling points the smaller the paramter_step
%      can be between iterations. Typical values to get started are 1000 and 0.05.
% (2)  This version of TRIXIFT can cope with scans of H, K, L or E, but not mixed scans.
%      In the event of scans involving multiple reciprocal lattice coordinates, edit 
%      the file 'trixfit.m'
% (3)  Before fitting EACH AND EVERY data set it is necessary to call trixfit_ini 
%
% Des McMorrow and Henrik Roennow
% Version: August 2002
%

initialise.monitor_flag=1;                   % 1 for count on monitor, 0 for count on time
initialise.monte_carlo_samples=2000;         % Monte Carlo steps for convolution integration    

initialise.resolution_method='rc_popma';     % rc_popma (Popovici) or rc_cnmat (Cooper-Nathans)
initialise.rescal_pars='lsco.par';           % parameters for Cooper-Nathans (mandatory)
initialise.popovici_pars='lsco.cfg';         % parameters for Popovici (optional)

initialise.xsec_file='lsco_xsec';            % definition of S(Q,w)
initialise.bkgd_file='lsco_bkgd';            % background definition 
initialise.pnam_file='lsco_pnam';            % sets parameter names
initialise.corr_file='lsco_corr';            % correction to calculated intensities, 
                                             % e.g. lambda/2 in monitor

%----- Run trifit_ini to initialise the parameters

error_status=trixfit_ini(initialise);
if ~isempty(error_status), disp('Error initialising parameters'),return, end

%----- Define initial value of parameters  

%----- Parameters specifying Qh,Qk,Ql and w for the scan, 
%      1000 indicates that it is the scan variable
                    
pin=[1000 0 -0.12 3];

%----- Cross-sec parameters defined in xsec_file
                    
delta=0.103;
pin=[pin 1.6 1 0 10 0.035 delta];

%----- Background parameters used in bkgd_file

pin=[pin 1 7]; 

%----- 3 meV data, incom. peak

% QH scan

%----- Initialise trixfit, load and combine data

  error_status=trixfit_ini(initialise);
  if ~isempty(error_status), disp('Error initialising parameters'),return, end
  s=loads('illbatch','data/0178[47 48 49],X=QH,Y=CNTS,M=M1')*135000;
  s=combine(0.005,s);

%----- Set fit control parameters (fcp), determine which parameters to fit and perform fit  
  
  fcp=[0.02 15 0.0001];
  notfixed=zeros(size(pin));
  notfixed(8:end)=1;
  [r,f]=fits(s,'trixfit',pin,notfixed,fcp);
  
%----- Plot data and fit  
  
  plot(r)  
  formatpars(f)
  xlabel('QH')
  ylabel('Counts per min=1.35e5')
  title('E=3 meV, kf=2.662 1/Å, T=1.6 K, 40min kf coll')
  legend(['Kappa=' num2str(pin(9))],2)
  axis([-inf inf 0 35])

%----- Plot starting guess  
  
  [x,y]=extract(s);
  f0=feval('trixfit',x,pin);
  line(x,f0,'linestyle','-.','color','r','linewidth',2)
  
  


