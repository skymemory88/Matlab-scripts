function [y, name, pnames, pin]=fitstr(x, p, flag)

hklc =[1     0     0
    -1     0     0
     0     1     0
     2     0     0
     0     2     0
    -2     0     0
     1     1     0
     1    -1     0
    -1     1     0
     2     1     0
     1     2     0
     2    -1     0
    -1     2     0
    -2     1     0
     3     0     0
     0     3     0
     3     1     0
    -1     3     0
    -3     1     0
     1    -3     0
     3    -1     0
     2     2     0
     2    -2     0
    -2     2     0
     3     2     0
     2     3     0
    -2     3     0
    -3     2     0
     2    -3     0
     3    -2     0];
 
hklc(2,:)=[];
 
hkla=[0     0     1
     0     0     2
     0     0     3
     0     0     4
     0     0     5
     0     0     6
     0     0     7
     1     0    -7
     1     0    -1
     1     0     0
     1     0     1
     1     0     2
     1     0     3
     1     0     4
     1     0     5
     1     0     6
     1     0     7
     2     0    -5
     2     0    -4
     2     0    -3
     2     0    -2
     2     0    -1
     2     0     0
     2     0     1
     2     0     2
     2     0     3
     2     0     4
     2     0     5
    -1     0    -2
    -1     0    -1
    -1     0     0
    -1     0     1
    -1     0     2
    -1     0     3
    -1     0     4
    -1     0     5
    -1     0     6
    -1     0     7
    -2     0    -3
    -2     0    -2
    -2     0    -1
    -2     0     0
    -2     0     1
    -2     0     2
    -2     0     3
    -2     0     4
    -2     0     5
    -2     0     6];
hkl=[hklc; hkla];
xhkl=x+29; % code: -28:0 c-orientierung 1:48 a-orientierung


if nargin==2
            %q=zeros(length(x)/3,3);
            %for ind=1:length(x)/3
            %    q(ind,1)=x(ind*3-2);
            %    q(ind,2)=x(ind*3-1);
            %    q(ind,3)=x(ind*3);
            %end
            
            q=hkl(xhkl,:);
             
           % y=intber(q,p')'; 
            %rotierte Domaene:
            faktor=p(1); %gewichtung nuclear reflexe
            parameter=p(2:13)'; 
          %  parameter(8)=p(9)*(1-p(13));
           % parameter(12)=p(9)*p(13);
            y=intber(q,[faktor,parameter])';

            if length(p)==25 %2.Domaene gewichtung frei
               parameter2=p(14:25)';
               y=(intens+intber(q,[faktor,parameter2])')/2;
            end
            if length(p)==14 %2.Domaene gedreht 
               gewichtung=p(14); %verhältnis domaene 1 zu Domaene 2 is 1:gewichtung
               qt=[q(:,2),-q(:,1),q(:,3)];
               y=(y+gewichtung*intber(qt,[faktor,parameter])')/(1+gewichtung);
            end
            
           % y=y/2+intber(qt,p')'/2;
           
        %if length(p)>13
        %    y(ind)=y(ind)+p(14)*intber(q,p(15:26)); %2.Domaene
        %end
else
	% ---------- Return initializtion data -----------------------
	y=[];
	name='Power law (y=0 x>x0)';
	pnames=str2mat('Amplitude','X.offset','Exponent','Background');
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
