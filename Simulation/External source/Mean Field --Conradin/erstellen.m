function erstellen()

temperatur=[0:0.1:2.4];
m=length(temperatur); 
for i=1:m
    punkte(i,:,1)=[0:0.1:6];
    punkte(i,:,2)=0;
    punkte(i,:,3)=0;
end
n=size(punkte,2);
Swerte=zeros(m,n,3);
Stemp=zeros(1,n,3);
putemp=zeros(1,n,3);
save('werte_ohnehyp5', 'temperatur', 'punkte', 'Swerte')

for iter=1:m
    ttemp=temperatur(iter);
    putemp(1,:,:)=punkte(iter,:,:);
    Stemp=gitter(putemp,ttemp);
    Swerte(iter,:,:)=Stemp(1,:,:);
    save('werte_ohnehyp5', 'Swerte','-append')
    iter
    
end
mesh(punkte(:,:,1),temperatur,Swerte(:,:,3));
figure
contour(temperatur,punkte(1,:,1),Swerte(:,:,3)')
ylabel('H/ Tesla')
xlabel('T/ Kelvin')
title(['\langle J_z \rangle MF LiHoF_4 without Hyperfeinww.'])

%mesh(punkte(:,:,1),temperatur,Swerte(:,:,3))

%temperatur=t;
%Swerte=S;
%punkte=pu;
 %Datei werte muss existieren '-append'
%clear temperatur Swerte punkte

%mesh(pu(:,:,1),t,S(:,:,1))
%mesh(pu(:,:,1),t,S(:,:,3))

if 1==1
pudx=punkte;
pudz=punkte;
t=temperatur;
for i=1:length(t)
    pudx(i,:,1)=punkte(i,:,1)+0.01;
    pudz(i,:,3)=0.01;
end

g=2

Sdx=gitter(pudx,t);
Sdz=gitter(pudz,t);

g=3

chixx=(Sdx(:,:,1)-Swerte(:,:,1))/0.01;
chizz=(Sdz(:,:,3)-Swerte(:,:,3))/0.01;

g=4
figure
mesh(punkte(:,:,1),t,chixx)
figure
plot(chixx(1,:))
%save('werte', 'chixx', 'chizz','-append')
end
