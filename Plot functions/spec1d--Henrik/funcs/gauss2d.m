function [y, name, pnames, pin]=gauss2d(x,p, flag)
% gauss     : Gaussian
% [y, {name, pnames, pin}]=gauss(x,p, {flag}) 
% MFIT Gaussian fitting function
% p = [ Amp Centrex Centrey  witdh1..width4 BackG]

% Author:  MZ <mzinkin@sghms.ac.uk>
% Description: Gaussian
global vebene1 vebene2

if nargin==2;

   Wmatrix=[p(4),p(5);p(5),p(7)];
   y=zeros(size(x));
   punktx=vebene1-p(2);
   punkty=vebene2-p(3);
   for n=1:length(x)
        y(n)=p(1)*exp(-0.5*[punktx(n) punkty(n)]*Wmatrix*[punktx(n) punkty(n)]')+p(8);%+(x-p(2))*p(9)
   end
  
   %plot(vy)
   
%    for n=1:N
%        if size(x,1)>=size(x,2)
%            z(n,:)=x([(n*N-N+1):n*N]);
%        else
%            z(:,n)=x([(n*N-N+1):n*N]);
%        end
%    end
%    %y=p(1)*exp(-0.5*(((x-p(2))/p(3)).^2))+p(4)+(x-p(2))*p(5);
%    z=
%    y=zeros(size(x));
%    for n=1:N
%        if size(x,1)>=size(x,2)
%            y((n-1)*N+1:n*N)=z(n,:);
%        else
%            y((n-1)*N+1:n*N)=z(:,n);
%        end
%    end
   
   
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
