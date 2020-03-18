function [y, name, pnames, pin]=gauss(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
   muB=0.05788*11.6;

   T=x;J=p(1);B=p(2);a=p(3);
   chi0=2*sinh(2*muB*B./T)./(exp(J./T)+1+exp(2*muB*B./T)+exp(-2*muB*B./T));
   chi=chi0./(1+p(6)*chi0);
   y=a*chi-p(4)./T-p(5);
   
else
	y=[];
	name='gaussian';
%	pnames=str2mat('Amplitude','Centre','Width','r','t','J''');
	pnames=str2mat('J','B','A','para','bck.','J''');
	if flag==1, pin=[0 1 1 1 1 1]; else pin = p; end
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
