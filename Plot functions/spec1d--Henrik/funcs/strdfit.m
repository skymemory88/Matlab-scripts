function [y, name, pnames, pin]=strdfit(x, p, flag)


if nargin==2
    nd1=[1:2:length(p)];
    nd2=[2:2:length(p)];
    pf1=p(nd1);%modellparameter fuer domain 1
    pf2=p(nd2);%modellparameter fuer domain 2
    %pf2(1)=p(2)*p(1); %p(1) skalierung gesamt, p(2) domain2:domain1
    pf2(1)=p(1);
    y=strfit(x,pf1)+strfit(x,pf2);
    y=y./sqrt(1+p(2)*y);

else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Power law (y=0 x>x0)';
	pnames=str2mat('skalierung','X.offset','Exponent','Background');
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

return

