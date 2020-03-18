function [y, name, pnames, pin]=gausss(x,p, flag)
% gausss     : Gaussian with sloping background
% [y, {name, pnames, pin}]=gausss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG Slope]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
    y=p(4)+p(1)*exp(-0.5*((x-p(2))/p(3)).^2)+p(5)*(x-p(2));
else
	y=[];
	name='gaussian';
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
