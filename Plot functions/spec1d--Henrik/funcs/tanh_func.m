function [y, name, pnames, pin]=gauss(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
   y=p(1)*tanh((x-p(2))/p(3));
   if length(p>=4)
      y=y+p(4);
   end
   if length(p>=5)
      y=y+p(5)*(x-p(2));
   end
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
