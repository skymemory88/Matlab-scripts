function [y, name, pnames, pin]=arrhenius_str(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Arhenius fitting function
% p = [ omega_0 EB ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Arrhenius Law

if nargin==2;
   y=p(1)-p(2)*x;
else
y=[];
    name='Arrhenius_str';
    pnames=str2mat('log(Omega zero)','Blocking Energy');
	if flag==1, pin=[1 1]; else pin = p; 
    end
end
