function [y, name, pnames, pin]=gauss_plus_voigt_amp(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ AmpGauss Centre WidthGauss AmpLorenz WidthLorenz BackGround Slope]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
   p(1:5)=abs(p(1:5));
   if length(p)<6
      p(6)=0;
   end
   if length(p)<7
      p(7)=0;
   end
   y=p(1)/(sqrt(2*pi)*p(3))*exp(-0.5*(((x-p(2))/p(3)).^2))+p(6)+(x-p(2))*p(7);
 %  if p(5)>0.5*p(3)
       lorzwidth=p(5)+0.01;
       if lorzwidth>10
           %lorzwidth=10;
       end
 %  else
 %      lorzwidth=0.5*p(3);
 %  end
  % if p(5)>0.5
  %     lorzwidth=0.5;
  % end

   y=y+voigt(x,[p(4) p(2) p(3) lorzwidth 0]);
   
       
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
