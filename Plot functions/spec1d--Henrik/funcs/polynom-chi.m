function [y, name, pnames, pin]=polynomchi(x,p, flag)
% polynomial: polynomial
% [y, {name, pnames, pin}]=polynomial(x,p, {flag})
%
% MFIT polynomial function 
%  y = p(1)*x^d + p(2)*x^(d-1) + ... + p(d)*x + p(d+1) p(5) = centre

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  polynomial
if nargin==2;
   y=p(1)*(x-p(5))^3 + p(2)*(x-p(5))^2 + p(3)*(x-p(5)) +p(4));
else
	y=[];
	name='gaussian';
   pnames=str2mat('Amplitude','Centre','Width');
   if length(p)>3
      pnames=str2mat(pnames,'Bck.');
   end
   if length(p)>4
      pnames=str2mat(pnames,'Slope');
   end
	if flag==1, pin=[0 1 1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width bg];
	end
end
