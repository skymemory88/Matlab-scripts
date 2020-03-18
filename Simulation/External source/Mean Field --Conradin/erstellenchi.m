function erstellenchi()

temperatur=[0:0.05:2];
m=length(temperatur); 
for i=1:m
    punkte(i,:,1)=[0:0.05:3.3];
    punkte(i,:,2)=0;
    punkte(i,:,3)=0;
end
n=size(punkte,2);
Swerte=zeros(m,n,3);
Stemp=zeros(1,n,3);
chiwerte=zeros(m,n,3);
chitemp=zeros(1,n,3);
putemp=zeros(1,n,3);
save('werte_ohnehyp4', 'temperatur', 'punkte', 'Swerte', 'chiwerte')

for iter=1:m
    ttemp=temperatur(iter);
    putemp(1,:,:)=punkte(iter,:,:);
    [Stemp,chitemp]=gitterchi(putemp,ttemp);
    Swerte(iter,:,:)=Stemp(1,:,:);
    chiwerte(iter,:,:)=chitemp(1,:,:);
    save('werte_ohnehyp4', 'Swerte','chiwerte','-append')
    iter   
end
mesh(punkte(:,:,1),temperatur,Swerte(:,:,3));
%figure
%contour(temperatur,punkte(1,:,1),Swerte(:,:,3)')
xlabel('H/ Tesla')
ylabel('T/ Kelvin')
title(['\langle J_z \rangle MF LiHoF_4 without Hyperfeinww.'])

figure
mesh(punkte(:,:,1),temperatur,chiwerte(:,:,3));
xlabel('H/ Tesla')
ylabel('T/ Kelvin')
title(['chi_zz MF LiHoF_4 without Hyperfeinww.'])


figure
mesh(punkte(:,:,1),temperatur,chiwerte(:,:,1));
xlabel('H/ Tesla')
ylabel('T/ Kelvin')
title(['chi_zz MF LiHoF_4 without Hyperfeinww.'])

%mesh(punkte(:,:,1),temperatur,Swerte(:,:,3))

%temperatur=t;
%Swerte=S;
%punkte=pu;
 %Datei werte muss existieren '-append'
%clear temperatur Swerte punkte

%mesh(pu(:,:,1),t,S(:,:,1))
%mesh(pu(:,:,1),t,S(:,:,3))

if 1==0
pudx=pu;
pudz=pu;

for i=1:length(t)
    pudx(i,:,1)=pu(i,:,1)+0.01;
    pudz(i,:,3)=0.01;
end

g=2

Sdx=gitter(pudx,t);
Sdz=gitter(pudz,t);

g=3

chixx=(Sdx(:,:,1)-S(:,:,1))/0.01;
chizz=(Sdz(:,:,3)-S(:,:,3))/0.01;

g=4

%mesh(pu(:,:,1),t,chixx)

save('werte', 'chixx', 'chizz','-append')
end
return
