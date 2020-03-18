function [y, name, pnames, pin]=polynomchi(x,p, flag)
% polynomial: polynomial
% [y, {name, pnames, pin}]=polynomial(x,p, {flag})
%
% MFIT polynomial function 
%  y = p(1)*x^d + p(2)*x^(d-1) + ... + p(d)*x + p(d+1), p(end) = centre

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  polynomial
if nargin==2;
    switch length(p)
        case 6
            y=p(1)*(x-p(end)).^4 + p(2)*(x-p(end)).^3 + p(3)*(x-p(end)).^2 +p(4)*(x-p(end)) +p(5);
        case 5 
            y=p(1)*(x-p(end)).^3 + p(2)*(x-p(end)).^2 + p(3)*(x-p(end))+p(4);
    end
else
	y=[];
	name='polynomial';
   pnames=str2mat('p1','p2','p3','p4','p5','centre');
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
