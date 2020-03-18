function [y, name, pnames, pin]=nvoigt(x,p, flag)
% ngauss    : N gaussians
% [y, {name, pnames, pin}]=nvoigt(x,p, {flag})
% 
% MFIT N voigt fitting function with common gaussian width
% p = [ Amp1 Centre1 Width1 ... AmpN CentreN WidthN Gauss_width Background [slope]]

% Author:  HMR <hmr@ill.fr>
% Description:  N voigts

ngss = ((length(p)-1)-mod((length(p)-1),3))/3;
if nargin==2
   
   if mod(length(p),3)<1
     p=[p(:);0];
   end
   if mod(length(p),3)<2
     p=[p(:);0];
   end

    y=p(end-1)+(x-p(2))*p(end);
    for i=1:ngss
        off=3*(i-1);
        y=y+voigt(x,[p(1+off) p(2+off) p(end-2) p(3+off) 0]);
    end;
else
   if flag==2
   	ngss=0;

      	mf_msg('Click on background');    
	[x bg]=ginput(1);

	but=1;
	pin=[];
		
	while but==1
		mf_msg(sprintf('Click on peak %d (right button to end)',ngss+1));
		[cen amp but]=ginput(1);
		if but==1
			mf_msg(sprintf('Click on width %d (right button to end)',ngss+1));
			[width y but]=ginput(1);
			width=abs(width-cen);
			amp=amp-bg;
			pin=[pin amp cen width];
			ngss=ngss+1;
		end
	end
	pin=[ pin bg];
   end

   y=[];
   pnames=str2mat('Amplitude_1','Centre_1','Width_1');
   if ngss>0
	name=sprintf('ngauss');
   else
	name='n gauss  : clik on Guess button to set n.';
	if flag==1, pin=[ 0 0 1 0 ]; else pin = p; end
   end

   for i=2:ngss
	pnames=str2mat(pnames,...
               sprintf('Amplitude_%d',i),...
               sprintf('Centre_%d',i),...
               sprintf('Width_%d ',i));
   end;
   if mod(length(p),3)>0
     pnames=str2mat(pnames,'Background');
   end
   if mod(length(p),3)>1
     pnames=str2mat(pnames,'Slope');
	end
  
end
