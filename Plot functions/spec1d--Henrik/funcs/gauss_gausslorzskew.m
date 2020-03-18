function [y, name, pnames, pin]=lorzskew(x,p, flag)
% lorz : Lorentzian
% [y, {name, pnames, pin}]=lorz(x,p, {flag})
%
% MFIT Skew Lorentzian fitting function
% p = [ Amplitude Centre Width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian

if nargin==2
   if length(p)<8
     p=[p(:);0];
   end
   if length(p)<9
     p=[p(:);0];
   end
    y=p(1)*exp(-0.5*(x-p(2)).^2/p(3)^2)+p(4)*(...
    (x<p(5)).*exp(-0.5*(x-p(5)).^2/p(6)^2)+...
    (x>=p(5))./(1+ (x-p(5)).^2/p(7)^2 ));
    y=y+p(8)+x*p(9);
else
	y=[];
	name='Lorentzian';
	pnames=str2mat('Amplitude','Centre','Width1','Width2','Background');
	if flag==1, pin=[0 0 1 1 1]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on width1');
		[width1 y]=ginput(1);
		width1=abs(width1-cen);
		mf_msg('Click on width2');
		[width2 y]=ginput(1);
		width2=abs(width2-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen width1 width2 bg];
	end
end
