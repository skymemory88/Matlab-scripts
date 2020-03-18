function [error_status]=trixfit_ini(initialise)
%
% TRIXFIT function to initialise all parameters
% 
% xsec:               cross-section file
% method:             rc_cnmat for Cooper-Nathans, rc_popma for Popovici
% parameter_file:     rescal parameter file
% configuration_file: configuration file for Popovici, can be specified as [] if 
%                     rc_cnmat is selected
%

global global_trixfit
global first_call

error_status=[];

global_trixfit.monitor_flag=initialise.monitor_flag;                   
global_trixfit.monte_carlo_samples=initialise.monte_carlo_samples;            

global_trixfit.bkgd_file=initialise.bkgd_file;
global_trixfit.pnam_file=initialise.pnam_file;
global_trixfit.corr_file=initialise.corr_file;

parameter_file=initialise.rescal_pars;
configuration_file=initialise.popovici_pars;
xsec=initialise.xsec_file;
method=initialise.resolution_method;

%----- Set defaults for input

if nargin==0
   disp('Must specify cross-section')   
   error_status=1;
   return
end

%----- Set initial value of flag used in triplefit

first_call=1;

%----- Initialise global parameters

pcn=[]; ppop=[];

%----- Read in CN parameters

fid=fopen(parameter_file);
if fid<0
   warning('Cooper-Nathans parameter file not in path')
   error_status=1;
   return
end
pcn=parameter_read(fid);
if length(pcn)~=42
   warning('Cooper-Nathans requires 42 parameters!')
   error_status=1;
   return
end

%----- Read in Popovici parameters

if isempty(configuration_file) & strcmp(method,'rc_popma')
   warning('Configuration file must be specified for Popovici method')
   error_status=1;
   return
end

if ~isempty(configuration_file) & ~strcmp(method,'rc_cnmat')
   fid=fopen(configuration_file);
   if fid<0
      warning('Popovici parameter file not in path')
      error_status=1;
      return
   end
   ppop=parameter_read(fid);
   if length(ppop)~=27
      warning('Popovici method requires 27 parameters!')
      error_status=1;
      return
   end
else
   ppop=zeros(27,1);   
end

%----- Define global variables

%----- Check to see if cross-section exists

check_method=exist(xsec);
if check_method~=2
   warning('Cross-section not in path')
   error_status=1;
   return
end

%---- Check to see if method exists

check_method=exist(method);
if check_method~=2
   warning('Resolution calculation method not in path')
   error_status=1;
   return
end

global_trixfit.xsec_file=xsec;
global_trixfit.resolution_method=method;
global_trixfit.pres=[pcn; ppop];
global_trixfit.first_call=first_call;

%=======================================================
function [p]=parameter_read(fid)
%
%
%

%-------------- Initialize arrays-----------------

data=[];		
header='';
text=fgetl(fid);

%----------------- Load data-----------------------

while (text>0)
   [temp count]=sscanf(text,'%f');			
   if isempty(temp) 
      header=[header text];
   else
      if (count==size(data,2) | isempty(data))
         data=[data; temp'];
      end
   end
      text=fgetl(fid);
end
fclose(fid);

p=data(:,1);
