function [y, name, pnames, pin]=gauss_plus_voigt(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ AmpGauss Centre WidthGauss AmpLorenz WidthLorenz BackGround Slope]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
   
   y=p(1)*exp(-0.5*((x-p(2))/p(3)).^2)+p(4);
 %  if p(5)>0.5*p(3)
       %lorzwidth=p(5)+0.01;
       %if lorzwidth>10
       %    lorzwidth=10;
       %end
 %  else
 %      lorzwidth=0.5*p(3);
 %  end
  % if p(5)>0.5
  %     lorzwidth=0.5;
  % end

   y=y+voigtn(x,[p(5) p(2) p(3) p(6) 0]);
   
       
else
	y=[];
	name='gaussian';
   pnames=str2mat('Amplitude','Centre','Width','bckg','area','lorzwidth');
   
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
