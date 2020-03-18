function [pnames]=trixfit_pnames(p)
%
% TRIXFIT function to define the fitting parameter names
%
% Des McMorrow, June 2001
%

pnames=strvcat('Qh','Qk','Ql','w');
pnames=strvcat(pnames,'T(K)','zc_H','zc_L','Int.','Zeta','Delta','Bckgd','Slope');

if length(pnames)~=length(p)
   warning('Number of parameter names not equal to number of parameters')
end