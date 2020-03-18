function [y, name, pnames, pin]=exponential(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Exponential fitting function
% p = [ Amp Centre Width BackG]

if nargin==2;
%    y=p(4)./(p(4)/p(1)*exp(-(x-p(2))/p(3))+1);
%    y=1./(exp(-(x-p(2))/p(3))/p(1)+p(4));
%    y=1./((p(2)-x).*exp((p(2)-x)/p(3))/p(1)+p(4));

    y=1./((p(2)-x).*exp((-x)/p(3))/p(1)+p(4));

else
	y=[];
	name='gaussian';
	pnames=str2mat('Amplitude','Centre','Width','Background');
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
