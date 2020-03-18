function [y, name, pnames, pin]=exponential(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Exponential fitting function
% p = [ Amp Centre Width BackG]

if nargin==2;
   p(1)=abs(p(1))/1000;
   p(2)=abs(p(2))/1000;
   x=x/1000;
   y=p(7)*x*1000-p(3)*(x<-p(1)).*(-x-p(1)).^p(5)+p(4)*(x>p(2)).*(x-p(2)).^p(6);
   if length(p)>7
      y=y+p(8);
   end
else
	y=[];
	name='gaussian';
	pnames=str2mat('V-','V+','c-','c+','n-','n+','d');
	if flag==1, pin=[0 0 1]; else pin = p; end
	if flag==2
		mf_msg('Sorry, no mouse input');
%		[cen amp]=ginput(1);
		pin=[1 1 1];
	end
end
