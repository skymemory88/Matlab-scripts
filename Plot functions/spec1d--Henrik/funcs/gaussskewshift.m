function [y, name, pnames, pin]=gaussskewshift(x,p, flag)
% gauss : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centre shift Width_left Width_r BackG ]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian

if nargin==2;
  sigma=p(4)*(x<=p(2))+p(5)*(x>p(2));
  y1=p(6)+p(1)*exp(-0.5*((x-p(2))./sigma).^2);
  sigma=p(4)*(x<=(p(2)+p(3)))+p(5)*(x>(p(2)+p(3)));
  y2=p(1)*exp(-0.5*((x-p(2)-p(3))./sigma).^2);
  y=y1-y2;  
else
	y=[];
	name='gaussian';
	pnames=str2mat('Amplitude','Centre','shift','Wl','Wr','Background');
	if flag==1, pin=[0 0 0 1 1 0]; else pin = p; end
	if flag==2
		mf_msg('Click on peak');
		[cen amp]=ginput(1);
		mf_msg('Click on shift');
		[shift y]=ginput(1);
      shift=shift-cen;
      mf_msg('Click on left width');
		[width y]=ginput(1);
		wl=abs(cen-width);
		mf_msg('Click on left width');
		[width y]=ginput(1);
		wr=abs(width-cen);
		mf_msg('Click on background');
		[x bg]=ginput(1);
		amp=amp-bg;
		pin=[amp cen shift wl wr bg];
	end
end
