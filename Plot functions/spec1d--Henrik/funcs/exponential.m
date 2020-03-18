function [y, name, pnames, pin]=exponential(x,p, flag)
% exponential     : Exponential
% [y, {name, pnames, pin}]=exponential(x,p, {flag}) 
% MFIT Exponential fitting function
% p = [ Amp Offset exponent BackG]

if nargin==2;
    y=p(4)+p(1)*exp((x-p(2))*p(3));
else
	y=[];
	name='exponential';
	pnames=str2mat('Amplitude','Offset','Exponent','Background');
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
