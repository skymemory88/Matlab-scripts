function [y, name, pnames, pin]=lorz(x,p, flag)
% lorz : Lorentzian
% [y, {name, pnames, pin}]=lorz(x,p, {flag})
%
% MFIT Lorentzian fitting function
% p = [ Amplitude Centre Width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian

if nargin==2
	y=p(4)+p(1)*p(3)./(p(3)^2 + (x-p(2)).^2 ) +p(5)*(x-p(2));
else
	y=[];
	name='Lorentzian';
	pnames=str2mat('Amplitude','Centre','Width','Background','Slope');
	if flag==1, pin=[0 0 1 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width');
		[width y]=ginput(1);
		width=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width bg 0];
	end
end
