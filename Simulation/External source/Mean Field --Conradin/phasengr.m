function phasengr()

load('werte_ohnehyp2', 'temperatur', 'punkte', 'Swerte')

smatrix=squeeze(Swerte(:,:,3))';
[tmatrix,hmatrix]=meshgrid(temperatur,punkte(1,:,1));

%mesh(tmatrix,hmatrix,smatrix) 
tmax=max(temperatur);
hmax=max(punkte(1,:,1));
intervall=[0:0.005:1];

for n=1:21
winkel=(n-1)*pi/40;
ta=tmax*intervall*cos(winkel);
ha=hmax*intervall*sin(winkel);
sa = diag(interp2(tmatrix,hmatrix,smatrix,ta,ha','cubic'));
xi=intervall(find(sa>0.05));
sai=sa(find(sa>0.05));
%yi = interp1(sai,xi,[0:0.01:4],'cubic');
%plot(intervall,sa,'.',xi,sai,'*',yi,[0:0.01:4])
x0=interp1(sai,xi,0.2,'cubic');
hwert(n)=hmax*sin(winkel)*x0;
twert(n)=tmax*cos(winkel)*x0;
end
plot(twert,hwert)
return
hw=zeros(length(temperatur));
for ind=1:length(temperatur)
    hw1=squeeze(punkte(ind,:,1));
    Sw1=squeeze(Swerte(ind,:,3));
    ih=find(Sw1<0.3);
  
    irp=1:ih(1)+10;
    if ih(1)>3
        ira=1:ih(1)-1;
        xi = interp1(Sw1(ira),hw1(ira),[0:0.01:1],'spline');
        figure
        plot(hw1([irp]),Sw1(irp),'.',xi, [0:0.01:1])
    else
        figure
        plot(hw1([irp]),Sw1(irp),'.')
    end
        
    %hw(ind)=hw1();
end
%plot(temperatur, hw)