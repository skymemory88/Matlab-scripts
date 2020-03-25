function Skw=qconst_wh(q,epsilon,t)
%berechnet chi(q,omega). q ist eine vektor, hvek entspricht Hx Werte, t constant.  

%load('werte_ohnehyp2', 'temperatur', 'punkte', 'Swerte')

%for ind=1:length(temperatur)
%    if temperatur(ind)==t
%        break
%    end
%end

%range=1:121;
%hvek=squeeze(punkte(ind,range,1));
%jxv=squeeze(Swerte(ind,range,1));
%jyv=squeeze(Swerte(ind,range,2));
%jzv=squeeze(Swerte(ind,range,3));
%clear temperatur punkte Swerte
hvek=[0:0.05:6];

N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N);
D(:,:,:,:)=dipole(q,15);
% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*(5/4)^2*0.05368;%*vol;

ow=[0.0001:0.005:0.5]; %omega
Hcf=cf();
chi0=zeros(3,3); %plus omeag
%chi0m=zeros(3,3); %minus omega

chi=zeros(3,3,length(ow),length(hvek));
%chim=zeros(3,3,length(ow),length(hvek));

for jh=1:length(hvek)
    [e,v]=mfieldz([0 hvek(jh) 0],t,1,0,1,Hcf,0);
    for jw=1:length(ow)
        chi0(:,:)=chi0_w([hvek(jh) 0 0],t,ow(jw),epsilon,e,v); 
        %chi0m(:,:)=chi0_w([hvek(jh) 0 0],t,-ow(jw),epsilon,e,v);
        M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
        Mm=zeros(3*N);
        for n=1:N  %n,m Summe �ber Ionen in der Einheitszelle
            for m=1:N
            M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:)*D(:,:,n,m);
         %  Mm((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0m(:,:)*D(:,:,n,m);
            end
        end
        chi(:,:,jw,jh)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
            ((eye(size(M))-M)\([chi0(:,:);chi0(:,:);chi0(:,:);chi0(:,:)]));
        %chim(:,:,jw,jh)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
        %    ((eye(size(Mm))-Mm)\([chi0m(:,:);chi0m(:,:);chi0m(:,:);chi0m(:,:)]));
    end
end


Skw=zeros(length(ow),length(hvek));
%Skwm=zeros(length(ow),length(hvek));
%Wirkungsquerschnitt berechnen
for nh=1:length(hvek)
    for nw=1:length(ow)
         if t==0
            Skw(nw,nh)=1/pi* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw,nh))));
         else
            Skw(nw,nh)= 1/pi/(1-exp(-ow(nw)*11.6/t))* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw,nh))));
               % Skwm(nw,nh)= 1/pi/(1-exp(-ow(nw)*11.6/t))* ...
               % sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                %    imag(chim(:,:,nw,nh))));
         end
    end 
end

save('disph_mithyp','hvek','ow','Skw');
figure
kleiner=(Skw<=1000);
groesser=1000*(Skw>1000);
pcolor(hvek,ow,log(Skw)); %.*kleiner+groesser
xlabel('H_x / Tesla')
ylabel('\omega / meV')
title(['T=', num2str(t), ',q =',mat2str(q),',\epsilon=0.01 - RPA LiHoF_4 without Hyperfeinww.'])
shading interp

mwechsel=[];
for nw=1:length(ow)
    maufab=0;
    for nh=2:length(hvek)
        if maufab==1
            if Skw(nw,nh)<=Skw(nw,nh-1)
                maufab=0;
                mwechsel=[mwechsel;[hvek(nh) ow(nw)]];
            end
        else
            if Skw(nw,nh)>=Skw(nw,nh-1)
                maufab=1;
            end
        end
    end
end
figure 
plot(mwechsel(:,1),mwechsel(:,2),'.')
xlabel('H_x / Tesla')
ylabel('\omega / meV')
title(['T=', num2str(t), ',q =',mat2str(q),',\epsilon=0.01 - RPA LiHoF_4 without Hyperfeinww.'])

%figure
%[maxw,Imax]=max(Skw');
%[maxw,Imax2]=max(Skw);
%plot(hvek(Imax),ow, hvek,ow(Imax2))

    