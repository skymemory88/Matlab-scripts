function S=gitter(punkt, t)
%punkt(i,j,k); i,j ist index für T bzw. H;
%dritter index k=1,2,3 nur falls H vektor
%Temperatur am punkt i,j ist t(i); HFeld am punkt i,j ist punkt(i,j,:)

%bsp für Argumente 
% t=[0:1:4]; for i=1:length(t), punkt(i,:)=[0:1:3]; , end


H=cf();
m=size(punkt,1); %ist immer gleich length(t)!
n=size(punkt,2); %Anzahl verschieden H-Punkte zur Temperatur t(i)

if ndims(punkt)==2  
    S=zeros(m,n,3); %falls nur Hx gegeben
    
    for i=1:m
        for j=1:n
            [sx,sy,sz]=remf(punkt(i,j),t(i),1,0,1,H);
            S(i,j,:)=[sx,sy,sz];
        end
    end
    
else
    S=zeros(m,n,3); %falls H vektor, d.h. punkte(i,j,:) entspricht [Hx,Hy,Hz] 
    
        for i=1:m
        for j=1:n
            hx=punkt(i,j,1);
            hy=punkt(i,j,2);
            hz=punkt(i,j,3);
            [sx,sy,sz]=remf([hx hy hz],t(i),1,0,1,H);
            S(i,j,:)=[sx,sy,sz];
        end
    end
end

%bsp Sx versus T ,H: surf(punkt,t,S(:,:,1)) oder falls H Vektor surf(punkt(:,:,1),t,S(:,:,1))
%save('werte', 'punkt', 'S', 't')
