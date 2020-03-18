function [y, name, pnames, pin]=glorz2(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
    y=gauss(x,[p(1),p(2),p(3),p(4),p(5)])+gauss(x,[p(1)*p(6),p(2)+p(7),p(3)*p(8),0])+lorz(x,[p(9),p(2),p(10),0]);
else
	y=[];
	name='gaussian';
   pnames=str2mat('Amplitude','Centre','Width');
   if length(p)>3
      pnames=str2mat(pnames,'Bck.');
   end
   if length(p)>4
      pnames=str2mat(pnames,'Slope');
   end
	if flag==1, pin=[0 1 1 1]; else pin = p; end
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
