function [x,y,err,xlab,ylab]=hmigui(filename)
%
%function [x,y,err,xlab,ylab]=hmigui(filename)
%
%   This function extracts scans from HMI data files
%   Ask user for 
%   the X,Y variables to import. 
%   EF 08.07.97 DFM 12.6.96 HMR 16.8.98
%
%Output variable
%   x        - the independent variable
%   y        - the dependent variable
%   err      - uncertainty in dependent variable = sqrt(y)
%   xlab     - name of 'x' data
%   ylab     - name of 'y' data

% uses : ffind.c as a mex file. illbatch.m

[x, y, err, xlab, ylab]=hmibatch([ filename ',gui' ]);




