function [y, name, pnames, pin]=lorzskew(x,p, flag)
% lorz : Lorentzian
% [y, {name, pnames, pin}]=lorz(x,p, {flag})
%
% MFIT Skew Lorentzian fitting function
% p = [ Amplitude Centre Width Background ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description:  Lorentzian

ngss = (length(p)-mod(length(p),4))/4;
if nargin==2
   if mod(length(p),4)<1
     p=[p(:);0];
   end
   if mod(length(p),4)<2
     p=[p(:);0];
   end

    y=zeros(size(x));
    for i=1:ngss
        off=4*(i-1);
        y=y+p(1+off)*(...
    (x< p(2+off)).*exp(-0.5*(x-p(2+off)).^2/p(3+off)^2)+...
    (x>=p(2+off)).*exp(-0.5*(x-p(2+off)).^2/p(4+off)^2));
    end;
    y=y+p(5+off)+(x-p(2+off))*p(6+off);
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
