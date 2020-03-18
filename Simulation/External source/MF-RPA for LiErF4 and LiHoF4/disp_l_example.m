% addpath('E:\2011\RPA_bastien')
% addpath('E:\2010\v2_LiErF4_march10\Spec1d\load')
% addpath('E:\2010\v2_LiErF4_march10\Spec1d\spec1d')
% clear all
%%%firsr domain
global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

temp=0.02;
omega=[-0.1:0.001:0.252];
epsilon=0.005;
h=0;
hvec=[0 0 h];
xEr=1;
exEr=0; % magnetic exchange in J*S_iS_j
exHo=-0.0001; %0.1 microeV from Henrik and Jens PRB
renorm_Ho=[1,1,1];
renorm_Er=[1,1,1];
withdemagn=true;
alpha=1; % shape is a sphere.
%Jex=0/((6/5)^2*0.05368);
%qh=[-1.2:0.1:-0.001 0.001:0.2];
% qh=[-1.2:0.1:-0.2];
% qvec=qh'*[0 0 1]+ones(size(qh'))*[0 0.002 0];
ql=[-1];
qvec=ql'*[0 0 1]+ones(size(ql'))*[0 0 0];
    % waitforbuttonpress()
% 
momente_er=[1 0 0
        -1 0 0
        -1 0 0
         1 0 0];
momente_ho=[0 0 1
        0 0 1
        0 0 1
         0 0 1];
% momente=[0 1 0
%          0 1 0
%          0 -1 0
%          0 -1 0];

%momente=remf(hvec,temp,momente,0,0.00043421);
momente=remf_2(hvec,temp,momente_er,momente_ho,exHo,exEr,xEr,1-xEr,0.23,0,true,alpha,renorm_Er,renorm_Ho);

%h,t,momente_er,momente_ho,exHo,exEr,xEr,xHo,ErHyp,HoHyp,withdemagn,alpha,renorm_Er,renorm_Ho
chi0=chi0_w(hvec,temp,momente,omega,epsilon,0,0.00043421);
%chi0m=chi0_w(hvec,temp,momente,-omega,epsilon,Jex);

for nq=1:size(qvec,1)
q=qvec(nq,:);
clear Skw
chi=chi_qw(q,chi0);
%chim=chi_qw(q,chi0m,Jex);

for nw=1:length(omega)
    if temp==0
        Skw(nw)=1/pi* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    else
        Skw(nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    end
end

srpa(nq)=spec1d(omega,Skw,sqrt(Skw))
%srpa(nq)=cut(srpa(nq),[-0.1 0.15]);
end

[sf(1),f(1)]=fits(cut(srpa(1),[-0.1 0.15]),'dho3',[1000,0.045,0.0052,73,0.13,1,0.022/11.6,0,73,0.089],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(1))
[sf(2),f(2)]=fits(cut(srpa(2),[-0.1 0.15]),'dho3',[1200,0.04,0.0052,50,0.125,1,0.022/11.6,0,95,0.1],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(2))
[sf(3),f(3)]=fits(cut(srpa(3),[-0.1 0.15]),'dho3',[1400,0.037,0.0052,200,0.13,1,0.022/11.6,0,40,0.12],[1 1 0 1 1 1 0 0 1 1]);
plot(sf(3))
[sf(4),f(4)]=fits(cut(srpa(4),[-0.1 0.15]),'dho3',[1500,0.037,0.0052,150,0.128,1,0.022/11.6,0,150,0.12],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(4))
[sf(5),f(5)]=fits(cut(srpa(5),[-0.1 0.15]),'dho3',[1500,0.04,0.0052,22,0.13,1,0.022/11.6,0,27,0.1],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(5))
[sf(6),f(6)]=fits(cut(srpa(6),[-0.1 0.15]),'dho3',[1400,0.045,0.0052,47,0.13,1,0.022/11.6,0,55,0.09],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(6))
[sf(7),f(7)]=fits(cut(srpa(7),[-0.1 0.15]),'dho3',[1200,0.05,0.0052,66,0.135,1,0.022/11.6,0,63,0.078],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(7))
[sf(8),f(8)]=fits(cut(srpa(8),[-0.1 0.15]),'dho3',[1000,0.058,0.0052,100,0.14,1,0.022/11.6,0,300,0.065],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(8))
[sf(9),f(9)]=fits(cut(srpa(9),[-0.1 0.15]),'dho3',[800,0.07,0.0052,130,0.14,1,0.022/11.6,0,300,0.06],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(9))
[sf(10),f(10)]=fits(cut(srpa(10),[-0.1 0.15]),'dho3',[900,0.08,0.0052,150,0.135,1,0.022/11.6,0,84,0.052],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(10))
[sf(11),f(11)]=fits(cut(srpa(11),[-0.1 0.15]),'dho3',[900,0.09,0.0052,120,0.13,1,0.022/11.6,0,33,0.045],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(11))
[sf(12),f(12)]=fits(cut(srpa(12),[-0.1 0.14]),'dho3',[800,0.1,0.0052,80,0.125,1,0.022/11.6,0,11,0.04],[1 1 0 1 1 1 0 0 1 1]);
plot(sf(12))
[sf(13),f(13)]=fits(cut(srpa(13),[-0.1 0.15]),'dho3',[1500,0.117,0.0052,150,0.12,1,0.022/11.6,0,150,0.038],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(13))
[sf(14),f(14)]=fits(cut(srpa(14),[-0.1 0.15]),'dho3',[1500,0.117,0.0052,150,0.12,1,0.022/11.6,0,5,0.036],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(14))
[sf(15),f(15)]=fits(cut(srpa(15),[-0.1 0.15]),'dho3',[800,0.1,0.0052,90,0.124,1,0.022/11.6,0,80,0.04],[1 1 0 1 1 0 0 0 1 1]);
plot(sf(15))


c1=[f(1).pvals(2) f(2).pvals(2) f(3).pvals(2) f(4).pvals(2) f(5).pvals(2) f(6).pvals(2) f(7).pvals(2) f(8).pvals(2) f(9).pvals(2) f(10).pvals(2)...
    f(11).pvals(2) f(12).pvals(2) f(13).pvals(2) f(14).pvals(2) f(15).pvals(2)]';
c2=[f(1).pvals(5) f(2).pvals(5) f(3).pvals(5) f(4).pvals(5) f(5).pvals(5) 0.13 f(7).pvals(5) f(8).pvals(5) f(9).pvals(5) f(10).pvals(5)...
    f(11).pvals(5) 0.124 f(13).pvals(5) f(14).pvals(5) f(15).pvals(5)]';
c3=[f(1).pvals(10) f(2).pvals(10) f(3).pvals(10) f(4).pvals(10) f(5).pvals(10) 0.09 f(7).pvals(10) f(8).pvals(10) f(9).pvals(10) f(10).pvals(10)...
    f(11).pvals(10) 0.04 f(13).pvals(10) f(14).pvals(10) f(15).pvals(10)]';
e1=spec1d(qvec(:,3),c1,0);
e2=spec1d(qvec(:,3),c2,0);
e3=spec1d(qvec(:,3),c3,0);
plot(e1,e2,e3)

framexy=[0.2 0.13 0.35 0.83];
a_handle=axes('position',framexy)
mapplot(srpa,ql,[],[])


[x,y,z]=mapplot(srpa,ql,[],[]);
  inds=find(-0.001<y<0.001);
 z(inds)=0;
  pcolor(x,y,z)

axis([-1.1 -0.1 -0.1 0.2])
 caxis([0,1000])
 shading flat
 xlabel('(0 -0.03 l) [r.l.u.]','fontname','helvetica','fontsize',16)
 ylabel('E (meV)','fontname','helvetica','fontsize',16)
 set(gca,'fontname','helvetica','fontsize',16,'xtick',[-1:0.5:-0.1])
 hold on
 plot(m1c,m2c)
 text(-0.9954,-0.0737,'RPA','fontname','arial','fontsize',14,'color','w')
 
%%%%%%%%%%%%second domain (result the same as the first domain)
global strategies; % global convergence strategies switches
strategies.powerlaw=true; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=true; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.

temp=0.022;
omega=[-0.1:0.001:-0.0001 0.0001:0.001:0.252];
epsilon=0.005;
h=0;
hvec=[0 0 h];
xEr=1;
exEr=0; % magnetic exchange in J*S_iS_j
exHo=-0.0001; %0.1 microeV from Henrik and Jens PRB
renorm_Ho=[1,1,1];
renorm_Er=[1,1,1];
withdemagn=true;
alpha=1; % shape is a sphere.
%Jex=0/((6/5)^2*0.05368);
%qh=[-1.2:0.1:-0.001 0.001:0.2];
% qh=[-1.2:0.1:-0.2];
% qvec=qh'*[0 0 1]+ones(size(qh'))*[0 0.002 0];
ql=[-1.25:0.1:0.15];
qvec=ql'*[0 0 1]+ones(size(ql'))*[0 -0.03 0];
    % waitforbuttonpress()
% 
momente_er=[0 1 0
        0 1 0
        0 -1 0
         0 -1 0];
momente_ho=[0 0 1
        0 0 1
        0 0 1
         0 0 1];

momente=remf_2(hvec,temp,momente_er,momente_ho,exHo,exEr,xEr,0,0.23,0,true,alpha,renorm_Er,renorm_Ho);

%h,t,momente_er,momente_ho,exHo,exEr,xEr,xHo,ErHyp,HoHyp,withdemagn,alpha,renorm_Er,renorm_Ho
chi0=chi0_w(hvec,temp,momente,omega,epsilon,0,0.00043421);

for nq=1:size(qvec,1)
q=qvec(nq,:);
clear Skw
chi=chi_qw(q,chi0);
%chim=chi_qw(q,chi0m,Jex);

for nw=1:length(omega)
    if temp==0
        Skw(nw)=1/pi* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    else
        Skw(nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    end
end

srpad(nq)=spec1d(omega,Skw,sqrt(Skw))
end
% 

for n=1:length(srpad)
    srpd(n)=cut(srpad(n),[-0.097 0.097]);
    srp2d(n)=cut(srpad(n),[0.105 0.2]);
    yyd=getfield(srpd(n),'y');
    xxd=getfield(srpd(n),'x');
    yy2d=getfield(srp2d(n),'y');
    xx2d=getfield(srp2d(n),'x');
    [mmd,indd]=max(yyd);
    [mm2d,ind2d]=max(yy2d);
    %gc(n)=cut(srp(n),[xx(ind)-0.02 xx(ind)+0.02]);
    %[mmg,indg]=max(yg);
    [sfd(n),fd(n)]=fits(srpad(n),'dho2',[mmd,xxd(indd),0.01,mm2d,xx2d(ind2d),2,0.022/11.6,0],[1 1 1 1 1 1 0 0]);
    %plot(srpa,srp,sf)
    %waitforbuttonpress()
   
end
for n=1:length(srpad)
    figure(n)
    plot(srpad(n),sfd(n))
end
ppd=[fd.pvals]';
%ppd(11,5)=0.136;
ppd(11,5)=0.135;

e1d=spec1d(qvec(:,3),ppd(:,2),0);
e2d=spec1d(qvec(:,3),ppd(:,5),0);
plot(e1d,e2d)

framexy=[0.2 0.13 0.35 0.83];
a_handle=axes('position',framexy)

mapplot(srpad,ql,[],[])
axis([-1.1 -0.1 -0.1 0.2])
 caxis([0,600])
 shading flat
 xlabel('(0 -0.03 l) [r.l.u.]','fontname','helvetica','fontsize',16)
 ylabel('E (meV)','fontname','helvetica','fontsize',16)
 set(gca,'fontname','helvetica','fontsize',16,'xtick',[-1:0.5:-0.1])
hold on
plot(m1c,m2c)


%%%%both domains
for m=1:length(srpad)
yj=getfield(srpa(m),'y');
xj=getfield(srpa(m),'x');
y2j=getfield(srpad(m),'y');
int=yj+y2j;
%et=e+e2;
stj(m)=spec1d(omega,int,0);
end

framexy=[0.2 0.13 0.35 0.83];
a_handle=axes('position',framexy)
axis([-1.1 -0.1 -0.1 0.2])
mapplot(stj,ql,[],[]);
caxis([0,1800])
shading flat
xlabel('(0 -0.03 l) [r.l.u.]','fontname','helvetica','fontsize',16)
 ylabel('E (meV)','fontname','helvetica','fontsize',16)
 set(gca,'fontname','helvetica','fontsize',16,'xtick',[-1:0.5:-0.1])
hold on
plot(m1c,m2c)

%%%%%%folded with exp res
global scalc_x scalc_y

for nn=1:length(srpa)
    scalc_x=getfield(srpa(nn),'x');
   
    scalc_y1=getfield(srpa(nn),'y');  scalc_y2=getfield(srpad(nn),'y');
    scalc_y=scalc_y1+scalc_y2;
%      yy=getfield(skn(nn),'y');
%      xx=getfield(skn(nn),'x');
%      [mm,ind]=max(yy);
     %gc(nn)=cut(st(nn),[x(ind)+0.01 0.3]);
%      yg=getfield(gc(nn),'y');
%      xg=getfield(gc(nn),'x');
%      [mmg,indg]=max(yg);
     
    %[sf(nn),f(nn)]=fits(sk(nn),'dho2_conv',[mm,x(ind),0.004,temp/11.6,0.02,0,0,mmg,xg(indg),0.012],[1 0 0 0 0 0 0 1 0 0]);
    %[sf,f]=fits(cut(sg(10),[-0.05 0.05]),'rpa_convolution',[0.01 1 1 1 0],[1 0 1 1 0]);%incoherentat 500mk y(inc)=1.7197
     [sl,fl]=fits(s(nn),'rpa_convolution',[0.017 1 0.0001 5.37 0],[1 0 1 1 0]);
     [sl,fl]=fits(s(nn),'rpa_convolution',[fl.pvals(1),1,1,fl.pvals(4:5)'],[0 0 1 0 1]);
     
    s_rpaconv(nn)=setfield(sl,'y',getfield(sl,'yfit'))
      
     %Plot(sf)
     plot(s_rpaconv(nn))
    % waitforbuttonpress()
end

framexy=[0.2 0.13 0.35 0.83];
a_handle=axes('position',framexy)
mapplot(s_rpaconv,ql,0.001,[],[]);
axis([-1.11 -0.1 -0.1 0.2])
 caxis([0,0.1])
 shading flat
 xlabel('(0 -0.03 l) [r.l.u.]','fontname','helvetica','fontsize',16)
 ylabel('E (meV)','fontname','helvetica','fontsize',16)
 set(gca,'fontname','helvetica','fontsize',16,'xtick',[-1:0.5:-0.1])
 hold on
 plot(m1c,m2c)
 axis([-1.11 -0.1 -0.1 0.2])
 

