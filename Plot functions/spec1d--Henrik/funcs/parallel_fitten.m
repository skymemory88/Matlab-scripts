function [y, name, pnames, pin]=parallel_fitten(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre Width BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

global fname xwerte separation pp sp ap

if nargin==2;
   y=zeros(size(x));
   for ns=1:length(separation)-1
       range=[separation(ns):separation(ns+1)-1];
       param=p(1+(ns-1)*sp:ns*sp);
       for m=1:length(pp)
           param(pp(m))=p(pp(m));
       end
       param(ap)=abs(param(ap));
       y(range)=feval(fname,xwerte(range),param);
   end
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
