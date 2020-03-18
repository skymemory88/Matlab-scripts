function [y, name, pnames, pin]=arrhenius(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Arhenius fitting function
% p = [ omega_0 EB ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Arrhenius Law

if nargin==2;
   y=p(1).*exp(-p(2)*x);
else
y=[];
    name='Arrhenius';
    pnames=str2mat('Omega zero','Blocking Energy');
	if flag==1, pin=[1 1]; else pin = p 
    end
end
