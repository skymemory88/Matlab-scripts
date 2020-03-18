function [y, name, pnames, pin]=lorz(x,p, flag)
% lorz      : Lorentzian
% [y, {name, pnames, pin}]=lorz(x,p, {flag})
%
% MFIT Lorentzian fitting function
% p = [ Amplitude Centre Width Background ]
% integral is : pi*p(1)*p(3) when bg = 0

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian

if nargin==2
	y=p(5)*p(6)./ (p(6)^2 + (x-p(2)).^2 )+p(1)*exp(-0.5*(((x-p(2))/p(3)).^2))+p(4);
else
	y=[];
	name='Lorentzian';
	pnames=str2mat('Gaussamp','Centre','Width','Background','lorzArea','lorzwidth');
	if flag==1, pin=[0 0 1 1]; else pin = p; end
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
