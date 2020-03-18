function [y, name, pnames, pin]=ngaussadd(x,p, flag)
% ngauss    : N gaussians
% [y, {name, pnames, pin}]=ngauss(x,p, {flag})
% 
% MFIT N gaussians fitting function
% p = [ Amp1 Centre1 Width1 ... AmpN CentreN WidthN Background]


% Author:  EF <manuf@ldv.univ-montp2.fr>, MZ <mzinkin@sghms.ac.uk>
% Description:  N gaussians

%zusaetzlich zu fittfunktion Korrektur
pp=p(4:end) %3 Korrekturparameter

%global ngss; 

ngss = (length(pp)-mod(length(pp),3))/3;
if nargin==2
   
   if mod(length(pp),3)<1
     pp=[pp(:);0];
   end
   if mod(length(pp),3)<2
     pp=[pp(:);0];
   end

    y=zeros(size(x));
    for i=1:ngss
        off=3*(i-1);
        y=y+pp(1+off)*exp(-0.5*(((x-pp(2+off))/pp(3+off)).^2));
    end;
    y=y+pp(4+off)+(x-pp(2+off))*pp(5+off);
    
    t=p(3);
    HoF3parameter =[0.7095   42.7825
    4.7600   16.8083
    4.0505    4.2117
   10.6646   11.0027];
    ef_HoF3=[2 3 3 4];
    ei_HoF3=[1 1 2 2];
    energie_Hof3 = [0
    0.7095
    4.7600
   11.3741
   11.8810
   20.8872
   21.9241
   23.5729
   25.9911
   28.7374
   29.2698
   31.5166
   32.8006
   33.9560
   41.7693
   44.5336
   45.1915];
    HoF3parameter(:,2) = p(1)*HoF3parameter(:,2).*exp(-energie_Hof3(ei_HoF3)/(t/11.6))/sum(exp(-energie_Hof3/(t/11.6)));
    HoF3parameter=[fliplr(HoF3parameter),p(2)*ones(size(HoF3parameter(:,1)))];
    HoF3parameter = reshape(HoF3parameter',1,numel(HoF3parameter));
    y=y+p(1)*ngauss(x,HoF3parameter)

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
