Gamma=0:0.01:1;
gamma=Gamma;
gamma(gamma<0.5)=gamma(gamma<0.5)*0+0.5;


wq1=sqrt(gamma.^2-0.5*1*Gamma.^2./gamma);
wq2=sqrt(gamma.^2-0.5*0.7*Gamma.^2./gamma);

q=[0:0.02:1];
Jq=cos(pi*q);

clear wq
for n=1:length(Jq)
  wq(n,:)=sqrt(gamma.^2-0.5*Jq(n)*Gamma.^2./gamma);
end

% Transverse Ising - dispersion
if 1==0
   clf
   set(gcf,'position',[100 100 420 280],'paperposition',[.25 2.5 6 4])
   axes('position',[0.15 0.15 .8 .8],'visible','off')
   axis([0 1 0 1])
   
   arrow([0 0],[1 0],20,'baseangle',60,'linewidth',2,'facecolor','k','edgecolor','k')
   arrow([0 0],[0 1],20,'baseangle',60,'linewidth',2,'facecolor','k','edgecolor','k')
   
  line(q,wq(:,(end+1)/2),'linewidth',2,'color','r')
  line(q,wq(:,1),'linewidth',2,'color','b')
  line(q,wq(:,round(5*(end+1)/8)),'linewidth',2,'color','b')

%  text(0.3,.15,'\tex[cc][cc]{$H=H_c$}','fontname','times','fontsize',2,'color','r')
%  text(0.1,0.6,'\tex[cc][cc]{$H=0$}','fontname','times','fontsize',2,'color','b')
%  text(0.6,0.9,'\tex[cc][cc]{$H>H_c$}','fontname','times','fontsize',2,'color','b')
%  text(0.9,-0.1,'\tex[cc][cc]{$q$}','fontname','times','fontsize',2)
%  text(-.15,0.95,'\tex[cc][cc]{$\omega_{q}$}','fontname','times','fontsize',2)
%  print -depsc -loose ising_disp.eps
  text(0.3,.15,'H=H_c','fontname','times','fontsize',12,'color','r')
  text(0.1,0.6,'H=0','fontname','times','fontsize',12,'color','b')
  text(0.6,0.9,'H>H_c','fontname','times','fontsize',12,'color','b')
  text(0.9,-0.1,'q','fontname','times','fontsize',12)
  text(-.15,0.95,'\omega_{q}','fontname','times','fontsize',12)
  print -depsc -loose ising_disp_2004.eps
end

% Transverse Ising - softening
if 1==0
   clf
   set(gcf,'position',[100 100 420 280],'paperposition',[.25 2.5 6 4])
   axes('position',[0.15 0.15 .8 .8],'visible','off')
   axis([0 1 0 1])
   arrow([0 0],[1 0],20,'baseangle',60,'linewidth',2,'facecolor','k','edgecolor','k')
   arrow([0 0],[0 1],20,'baseangle',60,'linewidth',2,'facecolor','k','edgecolor','k')
   
  line(Gamma,wq1,'linewidth',2,'color','r')
  line(Gamma,wq2,'linewidth',2,'color','b')

% text that gets converted by tex
%  text(0.7,0.15,'\tex[cc][cc]{$q=0$}','fontname','times','fontsize',3,'color','r','horizontalal','center')
%  text(0.4,0.55,'\tex[cc][cc]{$q\neq0$}','fontname','times','fontsize',3,'color','b','horizontalal','center')
%  text(0.5,-.1,'\tex[cc][cc]{$H_c$}','fontname','times','fontsize',3,'color','r','horizontalal','center')
%  text(0.9,-.1,'\tex[cc][cc]{$H$}','fontname','times','fontsize',3)
%  text(-0.15,.95,'\tex[cc][cc]{$\omega$}','fontname','times','fontsize',3)
%  print -depsc -loose ising_soften.eps
  text(0.7,0.15,'q=0','fontname','times','fontsize',12,'color','r','horizontalal','center')
  text(0.4,0.55,'q\neq0','fontname','times','fontsize',12,'color','b','horizontalal','center','interp','tex')
  text(0.5,-.1,'H_c','fontname','times','fontsize',12,'color','r','horizontalal','center','interp','tex')
  text(0.9,-.1,'H','fontname','times','fontsize',12)
  text(-0.15,.95,'\omega','fontname','times','fontsize',12,'interp','tex')
  print -depsc -loose ising_soften2004.eps
end


if ~exist('sz') & 1

T=0:0.01:0.5;
J=1;
gamma=zeros(size(T));
for n=1:length(T)
 if  T(n)==0
  gamma(n)=fzero_quiet(['x*',num2str(2/J,10),'-1'],[eps 10]);   
 elseif T(n)<0.25
  gamma(n)=fzero_quiet(['x*',num2str(2/J,10),'-tanh(x/2/',num2str(T(n),10),')'],[eps 10]);
 else
  gamma(n)=0;
 end
end
Gamma=0:0.025:1;

gamma=gamma'*ones(size(Gamma));
Gam=ones(size(T'))*Gamma;

sz=sqrt(max(gamma.^2-Gam.^2,0))/J;
sx=tanh(Gam./(T'*ones(size(Gamma)))/2)/2.*(sz==0);
tmp=Gam./gamma.*tanh(gamma./(T'*ones(size(Gamma)))/2)/2.*(sz>0);
sx(sz>0)=tmp(sz>0);
%sz(sqrt((T'*ones(size(Gamma))).^2+Gam.^2)>0.75)=sz(sqrt((T'*ones(size(Gamma))).^2+Gam.^2)>0.75)+nan;

end

if 1==1

clf
set(gcf,'position',[0 0 1120 420],'paperposition',[.25 2.5 16 6],'paperorientation','landscape')
axes('position',[0.05 0.1 0.3 0.7])
h1=mesh(T'*ones(size(Gamma)),Gam,sz);
view([70 30])
axis([0 .5 0 1 0 .5])
set(gca,'fontname','times','fontsize',24)
title('\langle S^z\rangle','fontsize',24,'fontname','times')
%xlabel('T/J','fontsize',24,'fontname','times')
%ylabel('\Gamma/J','fontsize',24,'fontname','times')
set(gca,'xticklabel',str2mat('0','T/J','1/2'),'xtick',[0 0.25 0.5])
set(gca,'yticklabel',str2mat('0','','1'))
text(0.65,0.5,'\Gamma/J','fontsize',24,'fontname','times','horizontalal','center')

axes('position',[0.35 0.1 0.3 0.7])
h2=mesh(T'*ones(size(Gamma)),Gam,sx);
view([70 30])
axis([0 .5 0 1 0 .5])
set(gca,'fontname','times','fontsize',24,'xticklabel','','yticklabel','','zticklabel','')
title('\langle S^x\rangle','fontsize',24,'fontname','times')

axes('position',[0.65 0.1 0.3 0.7])
h3=mesh(T'*ones(size(Gamma)),Gam,sqrt(sz.^2+sx.^2));
view([70 30])
axis([0 .5 0 1 0 .5])
set(gca,'fontname','times','fontsize',24,'xticklabel','','yticklabel','','zticklabel','')
title('|\langle{\bf{}S}\rangle|','fontsize',24,'fontname','times')

%print -depsc ising_mf_solution.eps
%!copy ising_mf_solution.eps c:\hmr\phd\lihof4

end