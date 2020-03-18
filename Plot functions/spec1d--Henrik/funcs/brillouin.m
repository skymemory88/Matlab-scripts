function [y, name, pnames, pin]=pow(x, p, flag)
% pow       : Power law y=0 x>x0
% function [y, name, pnames, pin]=pow(x,p, flag)
%
% MFIT Power law fitting function y=0 x>x0
% p = [ Amp Xoffset Exponent Background ]

% Author:  HR
% Description:  Brillouin function

if nargin==2
    S=p(3);
    Tc=p(2);
    M0=p(1);
    Bck=p(4);
    y=0*x;
    for n=1:length(x)
    T=x(n);
    z1=0;
    z2=100;
    for m=1:40
       z=(z1+z2)/2;
       B0=(1+1/(2*S))*coth((1+1/(2*S))*z)-1/(2*S)*coth(z/(2*S))-(S+1)/(3*S)*T/Tc*z;
       if B0<0
           z2=z;
       else
           z1=z;
       end
    end
    y(n)=(z1+z2)/2*(S+1)/(3*S)*T/Tc*M0+Bck;
end

else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Power law (y=0 x>x0)';
	pnames=str2mat('Mo','Tc','S','Bck');
	if flag==1, pin=[0 0 0 0]; else pin = p; end
   if flag==2
        %-------------- Get starting guesses 
        % Get starting guess for x offset
        mf_msg('Click on x offset estimate (y=0 x>x0)');	    
        [xoff dummy]=ginput(1);

        % Get starting guess for background 
        mf_msg('Click on background estimate');	    
        [dummy bg]=ginput(1);

        % Get two points on the 'curve'
        mf_msg('Click on two points on curve');	    
        [xp yp]=ginput(2);			   

        % Work out starting guess for the exponent
        power=log((yp(1)-bg)/(yp(2)-bg))/log((xp(1)-xoff)/(xp(2)-xoff));

        % Work out starting guess for the amplitude
        amp=abs((yp(1)-bg)/((xp(1)-xoff)^power));
    
        pin=[amp xoff power bg];
    end
end

