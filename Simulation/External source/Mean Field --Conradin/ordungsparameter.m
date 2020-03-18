h=[0:0.02:4];
t=[0.1:0.1:1];

Jx=zeros(length(h),length(t));
Jy=zeros(length(h),length(t));
Jz=zeros(length(h),length(t));
for n=1:length(h)
    for m=1:length(t)
        [Jx(n,m),Jy(n,m),Jz(n,m)]=remf(h(n),t(m),1,0,1);
    end
end

plot(h,Jz(:,1),h,Jz(:,2))

for n=1:length(t)
    [minimum,ind]=min(diff(Jz(:,n)));
    Hcr(n)=h(ind+1);
end
plot(t,Hcr)

[jx,jy,jz]=remf(2.9,0.2,1,0,1);
epsilon=0.01;
omega=[-0.3:0.01:1];
chi0=chi0_w([2.9 0 0],0.2,[jx,jy,jz],omega,epsilon);
q=[1,0,1];
t=0.2;

qv=[];
for l=0.5:0.05:1.5
qvt=[0:0.05:1]'*[1 0 0];
qvt(:,3)=l;
qv=[qv;qvt];
end
for nq=1:size(qv,1)
    q=qv(nq,:);
    chi=chi_qw(q,chi0);
for nw=1:length(omega)
if t==0
    Skw(nq,nw)=1/pi* ...
        sum(sum((eye(3)-(q'*q)/(q*q')).* ...
            imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
else
    Skw(nq,nw)= 1/pi/(1-exp(-omega(nw)*11.6/t))* ...
        sum(sum((eye(3)-(q'*q)/(q*q')).* ...
            imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
end
end
end

pcolor(Skw)


%%%
t=[0.1:0.1:1];
for n=1:length(Hcr)
    n
[jx,jy,jz]=remf(Hcr(n),t(n),1,0,1);
epsilon=0.01;
omega=[-0.3001:0.01:1];
chi0=chi0_w([Hcr(n) 0 0],t(n),[jx,jy,jz],omega,epsilon);
qv=[1.1,0,1
    1.05,0,0];
for nq=1:size(qv,1)
    q=qv(nq,:);
    chi=chi_qw(q,chi0);
for nw=1:length(omega)
if t(n)==0
    Skw(n,nq,nw)=1/pi* ...
        sum(sum((eye(3)-(q'*q)/(q*q')).* ...
            imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
else
    Skw(n,nq,nw)= 1/pi/(1-exp(-omega(nw)*11.6/t(n)))* ...
        sum(sum((eye(3)-(q'*q)/(q*q')).* ...
            imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
end
end
end
end
plot(omega,squeeze(Skw(1:5,2,:))')


