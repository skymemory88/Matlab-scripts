function [y, name, pnames, pin]=lorzn(x,p, flag)
% lorzn     : Lorentzian Power n
% [y, {name, pnames, pin}]=lorzn(x,p, {flag})
%
% MFIT Lorentzian power n fitting function
% p = [ Amplitude Centre Width Power Background ]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Lorentzian Power n

if nargin==2
%   y=p(5)+p(1)./abs((1+ (x-p(2)).^2/p(3)^2 ).^(p(4)+p(6)*x));
   y=p(5)+p(1)./abs((1+ (x-p(2)).^2/p(3)^2 ).^p(4))+...
          p(6)./abs((1+ (x-p(7)).^2/p(8)^2 ).^p(9));
   gN=real(1./(4*pi^2*sqrt(p(10)^2-x.^2+0.05i)));
   y=real(y.*gN);
   if p(11)>0
      y=y./abs(exp(-x/p(11))-1);
   else
      y=y.*(x>=0);
   end
else
	y=[];
	name='Lorentzian Power N';
%	pnames=str2mat('Amplitude','Centre','Width','Power','Background','Power2');
	pnames=str2mat('Amplitude','Centre','Width','Power','Background','Amplitude2','Centre2','Width2','Power2','q','Temp');
	if flag==1, pin=[0 0 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width 1 bg];
	end
end

