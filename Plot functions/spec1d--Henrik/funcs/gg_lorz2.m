function [y, name, pnames, pin]=gg_lorz(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2
    p(7)=abs(p(7));
    if abs(p(6)) > 1e4 || abs(p(6)) < 1e-4
        p(6) = 100;
    end
    if p(7) > 100 || p(7) < 1e-3
        p(7) = 2;
    end
    y=gauss(x,[p(1),p(2),p(3),p(4),p(5)])+lorz(x,[p(6),p(2),p(7),0])+gauss(x,[p(1)*p(8),p(9),p(3),0,0]);
else
	y=[];
	name='2gaussian+Lorentzian';
   pnames=str2mat('Amplitude','Centre','Width','Background','Slope','Amp','Width','gauss scaling','centre2');
	if flag==1, pin=[0 1 1 1 1 1 1]; else pin = p; end
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