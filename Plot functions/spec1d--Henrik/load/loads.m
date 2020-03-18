function s=loads(filetype,datafile)
%
% function s=loads(filetype,datafile)
%
% MATLAB function to load datafile of type filetype  
% to a spec1d object s.
%
% Filetype must be a valid load routine
% These include:  
%
% xyeload     :  x,y,error
% multibatch  :  multicolumn
% specbacth   :  SPEC file
% tasbatch    :  TASCOM file
% illbatch    :  ILL 3-axis file
% parbatch    :  MFIT parameter output file
% 
% Eaxmples  
%
% Load the 3rd scan  from SPEC file cesb.05
% >s=loads('specbatch','cesb.05,X=K,Y=Detector,M=Monitor,S=3');    
% 
% Load the TASCOM file cuge108.dat
% >s=loads('tasbatch','cuge108.dat,X=OM,Y=I,M=Mon');
%
% DFM 1.4.98
s=[];

[x,y,err,xlab,ylab]=feval(filetype,datafile);

if isempty(x); return; end

datafile=strtok(datafile,',');

s.x=x;
s.y=y;
s.e=err;
s.x_label=xlab;
s.y_label=ylab;
s.datafile=datafile;
s.yfit=[];

s=spec1d(s);
