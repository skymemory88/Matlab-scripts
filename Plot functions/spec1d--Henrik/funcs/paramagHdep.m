function [y, name, pnames, pin]=paramagHdep(x,p, flag)
% Author:  LQM
% Description: Paramagnetic susceptibility as fct. of field
% p(1) = amplitude
% p(2) = temperature 
% p(3) = background

if nargin==2;
   if length(p)<3
      p(3)=0;
   end
   y=p(1)*(1-tanh((x)./p(2)).^2)./p(2)+p(3);
else
	y=[];
	name='Paramag Hdep';
   pnames=str2mat('Amplitude','Temp');
   if length(p)>2
      pnames=str2mat(pnames,'Bck.');
   end
	if flag==1, pin=[0 1 1]; else pin = p; end
	if flag==2
		pin=[0 0 0];
	end
end
