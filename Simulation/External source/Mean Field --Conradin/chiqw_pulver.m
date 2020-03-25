function [fpar,sp,Nq]=chiqw_pulver(qbetrag,epsilon,t)



N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N);

Nqpunkte=100;


a=[5.175 0 0
   0 5.175 0
   0 0 10.75];
% Positions of the moments within the unit cell
tau=[0 0 0
     0 1/2 1/4
     1/2 1/2 1/2
     1/2 0 3/4];
 
 tau=tau*a;

% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];
% Convert q from reciprocal aangstroms to Miller indicies

fpar=[];
%nphi=[0;0];
%ntheta=[-pi; pi];
%delta=pi/10;
%for theta=(-pi/2+delta):delta:(pi/2-delta)
%    Nph=round(abs(cos(theta)/0.02));
%    for n=1:Nph
%        phi=n*pi*2/Nph;
%        nphi=[nphi;phi];
%        ntheta=[ntheta;theta];
%    end
%end
        

xyz=1.1*qbetrag*rand(1000,3);
betrg=sqrt(sum(xyz.^2,2));
qxyz=xyz(find(betrg<1.1*qbetrag&betrg>0.9*qbetrag),:);
Nq=size(qxyz(:,1));

%Nq=length(nphi);
for nwinkel=1:Nq
   
    
  %    [qx,qy,qz]=sph2cart(nphi(nwinkel),ntheta(nwinkel),qbetrag);
  %    q=[qx, qy, qz];
  q=[qxyz(nwinkel,1),qxyz(nwinkel,2),qxyz(nwinkel,3)];
      
D(:,:,:,:)=dipole(q*inv(b),15);
% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*(5/4)^2*0.05368;%*vol;

ow=[0.75:0.01:1.1]; %omega
Hcf=cf();
chi0=zeros(3,3); %plus omega
%chi0m=zeros(3,3); %minus omega

chi=zeros(3,3,length(ow));
%chim=zeros(3,3,length(ow),length(hvek));


    [e,v]=mfieldz(0,t,0,0,0,Hcf,1);
    for jw=1:length(ow)
        chi0(:,:)=CHI0_W(0,t,ow(jw),epsilon,e,v); 
        %chi0m(:,:)=chi0_w([hvek(jh) 0 0],t,-ow(jw),epsilon,e,v);
        M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
        Mm=zeros(3*N);
        for n=1:N  %n,m Summe �ber Ionen in der Einheitszelle
            for m=1:N
            M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:)*D(:,:,n,m);
         %  Mm((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0m(:,:)*D(:,:,n,m);
            end
        end
        chi(:,:,jw)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
            ((eye(size(M))-M)\([chi0(:,:);chi0(:,:);chi0(:,:);chi0(:,:)]));
        %chim(:,:,jw,jh)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
        %    ((eye(size(Mm))-Mm)\([chi0m(:,:);chi0m(:,:);chi0m(:,:);chi0m(:,:)]));
    end



Skw=zeros(length(ow),1);

%Wirkungsquerschnitt berechnen

    for nw=1:length(ow)
         if t==0
            Skw(nw)=1/pi* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));
         else
            Skw(nw)= 1/pi/(1-exp(-ow(nw)*11.6/t))* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw))));
               % Skwm(nw,nh)= 1/pi/(1-exp(-ow(nw)*11.6/t))* ...
               % sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                %    imag(chim(:,:,nw,nh))));
         end
    end 

sp(nwinkel)=spec1d(ow,Skw,sqrt(Skw));
[spft,ft]=fits(sp(nwinkel), 'lorz',[1,0.9,epsilon/2,0],[1 1 1 1]);
%plot(sp)
%waitforbuttonpress
%line([0:0.01:1.5],lorz([0:0.01:1.5],ft.pvals),'Color','r','LineWidth',1)
%line([2.5:0.01:4.2],lorz([2.5:0.01:4.2],[ft.pvals(1),ft.pvals(2),epsilon,ft.pvals(4)]) )
%waitforbuttonpress
fpar=[fpar;[ft.pvals(1:3)]'];
  
end





return
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

    