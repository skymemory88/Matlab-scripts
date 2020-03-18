function [y, name, pnames, pin]=dho(x,p, flag)
% dho       : Damped harmonic oscillator
% [y, {name, pnames, pin}]=dho(x,p, {flag})
%	
% MFIT  Damped harmonic oscillator.
% The function contains a Bose correction for the temperature factor
%  Physica B 234-236 (1997) 1107
%
% Units: T in 'x' unit for bose factor
% 1 meV = 242 GHz = 8.065 cm-1 = 11.604 K = 12.398 Angs
% DFM 19.10.95 rev EF 27.06.97
% p = [ DHO_Amp DHO_centre DHO_Half_Width T BackG ]

% Author:  EF <manuf@ldv.univ-montp2.fr>
% Description:  Damped harmonic oscillator

if nargin==2

    omega=x(:);
    amp=p(1);
    omega_0=p(2); 
    gamma=p(3);	
    T=p(4);

    n_omega=1./(1-exp(-omega/T));	% bose
    omega2 = omega.^2;


    y=4*amp*gamma*omega.*n_omega./((omega2-omega_0^2).^2+4*gamma^2*omega2)/pi+p(5);

else

	y=[];
	if flag==1, pin=[ 0 0 1 1 1]; else pin = p; end
	name='Damped harmonic oscillator';
	pnames=str2mat('DHO.Amplitude','DHO.Centre','DHO.Gamma',...
		       'T ("x" unit)','Background');

   if flag==2

      		mf_msg('Click on background');    
		[x bg]=ginput(1);

		mf_msg('Click on DHO Centre');
		[cen amp]=ginput(1);
		mf_msg('Click on width of DHO');
		[width y but]=ginput(1);
		width=abs(width-cen);
		amp=(amp-bg);
		cen=cen;
		pin=[amp/cen/cen cen width 1 bg ];

	end

end
