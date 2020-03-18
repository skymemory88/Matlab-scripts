function [y, name, pnames, pin]=powderpeaks(x,p, flag)
% powderpeaks
% [y, {name, pnames, pin}]=powderpeaks(x,p, {flag}) 
% MFIT fitting function
% x= list of tau's for the powder
% y= list of corresponding 2Theta values

% p = [wavevector 2Theta offset]

% HMR 2001

if nargin==2;
    y=2*asin(x/2/p(1))*180/pi-p(2);
else
	y=[];
	name='powderpeaks';
	pnames=str2mat('Wavevector','2Theta offset');
	if flag==1, pin=[1 0]; else pin = p; end
	if flag==2
		mf_msg('Sorry, no mouse input - GUI freak!');
		pin=[1 0];
	end
end
