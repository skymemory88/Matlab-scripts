function chiqw_schnitt(t)

%load('werte_ohnehyp2', 'Swerte', 'punkte', 'temperatur')

%for ind=1:length(temperatur)
%    if temperatur(ind)==t
%        break
%    end
%end
%szut=squeeze(Swerte(ind,:,3));
%igrenze=find(szut<0.1);
%h=punkte(ind,igrenze(1),:);
%jxyz =squeeze(Swerte(ind,igrenze(1),:));
%[e v]=mfieldz(h,t,jxyz(1),jxyz(2),jxyz(3),cf(),0);
h=0;%3.17;
[e v]=mfieldz(h,t,1,0,1,cf(),0);

qv=[1:0.01:2]'*[1 0 0];
ow=[0.0001:0.005:0.6];

%Skw=zeros(3,3,length(ow),size(qv,1));
Skw=zeros(length(ow),size(qv,1));


chi=chi_qw(qv,ow,0.01,h,t,e,v);
chim=chi_qw(qv,-ow,0.01,h,t,e,v);
for nq=1:size(qv,1)
    for nw=1:length(ow)
        q=qv(nq,:);
        if t==0
            Skw(nw,nq)=1/pi* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
        else
            Skw(nw,nq)= 1/pi/(1-exp(-ow(nw)*11.6/t))* ...
                sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                    imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
        end
    end
    qachse(nq)=qv(nq,1);
end

kleiner=(Skw<=4000);
groesser=4000*(Skw>4000);
figure
pcolor(qachse,ow,log(Skw))%.*kleiner+groesser
xlabel('[q,0,0]')
ylabel('\omega / meV')
title(['T=', num2str(t), ', H_x =',num2str(h(1)),',\epsilon=0.01 - RPA LiHoF_4 with Hyperfeinww.'])
shading flat


    