function phasengrenze()
   
%Phasengrenze, chixx maximal 

feldrichtung=1;


if feldrichtung==1
    strfr='x';
end


tmax=2; 
hmax=3.3;

intervall=[0.8:0.02:1];
intervall2=[-0.02:0.002:0.02];
winkel=[0:0.025:0.5]*pi;
winkel(length(winkel))=winkel(length(winkel))-0.001;
xy=zeros(length(intervall),2);
xy2=zeros(length(intervall2),2);
chixx=0*intervall;
chixx2=0*intervall2;
th=zeros(length(winkel),2);
H=cf();

for n=1:length(winkel)
    for m=1:length(intervall)
        [xy(m,1),xy(m,2)]=pol2cart(winkel(n),intervall(m));
        xy(m,:)=[tmax,hmax].*xy(m,:);
        hfeld=[0 0 0];
        hfeld(feldrichtung)=xy(m,2);
        [sx,sy,sz,chi_alphabeta]=remfchi(hfeld,xy(m,1),1,0,1,H);
        chixx(m)=real(chi_alphabeta(3,3)); 
    end
    [mw,mi]=max(chixx);
    intervall3=intervall(mi)+intervall2;
    for m=1:length(intervall3)
        [xy2(m,1),xy2(m,2)]=pol2cart(winkel(n),intervall3(m));
        xy2(m,:)=[tmax,hmax].*xy2(m,:);
        hfeld=[0 0 0];
        hfeld(feldrichtung)=xy2(m,2);
        [sx,sy,sz,chi_alphabeta]=remfchi(hfeld,xy2(m,1),1,0,1,H);
        chixx2(m)=real(chi_alphabeta(3,3)); 
    end
    [mw,mi]=max(chixx2);
    plot(intervall,chixx,'+',intervall3,chixx2,'.')
    figure
    th(n,:)=xy2(mi,:);
end

plot(th(:,1),th(:,2));
title(['Phasengrenze RE-Meanfield Antiferromagent LiHoF_4, Feld entlang ',strfr])
ylabel(['H_',strfr,' / Tesla'])
xlabel('T / Kelvin')


return
