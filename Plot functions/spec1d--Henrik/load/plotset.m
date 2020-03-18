if exist('lsmo03feb_220K_100AA_array.mat')
   load lsmo03feb_220K_100AA_array
   V=a.V;I=a.I;Vz=a.Vz;Iz=a.Iz;T=a.T;s=a.s;dIdV0=a.dIdV0;par=a.par;Is=a.Is
else
   fileroot='c:\STMdata\03mar13\13mar';
   files=[1:992 994:1985]; 
   [V,I,Vz,Iz,T,s]=stmdatas(fileroot,files);
   % Fit 5th order polynomial, slope is coefficient of linear term.
   [stmp,fstr]=fits(cut(s,[-300 300])*[0.001 1],'polynomial',[0 0 0 0 1 -0.014],[1 0 1 1 1 1]); % determine zero bias slope, and fix in global fits.
   par=[[fstr.pvals]' [fstr.evals]'];
   dIdV0=[par(:,5) par(:,11)];
   Is=smooth(V,I,Vsmooth);
   a.V=V;a.I=I;a.Vz=Vz;a.Iz=Iz,a.T=T;a.s=s;a.dIdV0=dIdV0;a.par=par;a.Is=Is;
   save lsmo03feb_220K_100AA_array a
end

% Fit Gaussian to Iz distribution
[Nh,Ih]=hist(Iz,0:0.005:1);
sum(Iz<0.3)/length(Iz)
sh=spec1d(Ih,Nh,sqrt(Nh));
[stmp,fh]=fits(cut(sh,[0 0.2]),'gauss',[250 .18 .02 0],[1 1 1 0]);

hist(Iz,0:0.005:1);axis([0 1 0 200])
nicefig
xlabel('I(Vz)')
ylabel('Promille')
title('LSMO 2003/3, T=220 K, Vz=0.8 V, Iz=0.2 nA')
line(Ih,gauss(Ih,fh.pvals),'linewidth',1,'color','r')
text(.4,175,sprintf('I_z= %3.3f \\pm %3.3f',[fh.pvals(2:3)]),'fontname','times','fontsize',16)
text(.4,150,sprintf('Total %4.4i',length(Iz)),'fontname','times','fontsize',16)
text(.4,125,sprintf('In Gaussian   %2.2i %%',round(sqrt(2*pi)*fh.pvals(3)*fh.pvals(1)/0.005/length(Iz)*100)),'fontname','times','fontsize',16)
text(.4,100,sprintf('Below I_z=0.3   %2.2i %%',round(sum(Iz<0.3)/length(Iz)*100)),'fontname','times','fontsize',16)
line([0.3 0.3],[0 50],'linewidth',1,'color','r','linestyle','--')
print -depsc lsmo03feb_220K_100AA_Izhist

% Try looking only at those spectra having I(Vz) close to the average
ngood=find(Iz>0.016 & Iz<0.2);
ngood=find(Iz>0);

[Nh,Vh]=hist(dIdV0(:,1),-0.003:0.000075:0.02);
sh=spec1d(Vh,Nh,sqrt(Nh));
[stmp,fh]=fits(cut(sh,[-1 0.2]),'gauss',[250 0 .001 0],[1 1 1 0]);


clf
hist(dIdV0(ngood,1),-0.003:0.000075:0.003)
line(Vh,gauss(Vh,fh.pvals),'linewidth',1,'color','r')
nicefig
axis([-0.0015 0.0025 0 120])
xlabel('dI/dV at V=0 [nA/V]')
ylabel('Promille')
title('LSMO 2003/3, T=220 K, Vz=0.8 V, Iz=0.2 nA, (100 Å)^2')
text(-1.4e-3,110,sprintf('dI/dV= %2.2f \\pm %2.2f pA',[fh.pvals(2:3)']*1000),'fontname','times','fontsize',12)
print -depsc lsmo03feb_220K_100AA_dIdV0hist

clf
hist(dIdV0(ngood,1),-0.003:0.00015:0.02)
line(Vh,gauss(Vh,fh.pvals)*2,'linewidth',1,'color','r')
nicefig
axis([-0.0015 0.02 0 240])
xlabel('dI/dV at V=0 [nA/V]')
ylabel('Promille')
title('LSMO 2003/3, T=220 K, Vz=0.8 V, Iz=0.2 nA, (100 Å)^2')
text(5e-3,210,sprintf('dI/dV= %2.2f \\pm %2.2f pA',[fh.pvals(2:3)']*1000),'fontname','times','fontsize',12)
text(5e-3,170,sprintf('In Gaussian   %2.2i %%',round(sqrt(2*pi)*fh.pvals(3)*fh.pvals(1)/0.000075/length(Iz)*100)),'fontname','times','fontsize',12)
print -depsc lsmo03feb_220K_100AAb_dIdV0hist

% Extract Gap
th=3e-4;
clear Vm Vp
for n=1:size(V,2)
   ng=find(Is(:,n)-V(:,n)/1000*par(n,5)-par(n,6) < -th);
   if isempty(ng)
      ng=1;
   end
   Vm(n)=-V(min(ng),n);
   Vp(n)= V(max(find(Is(:,n)-V(:,n)/1000*par(n,5)-par(n,6) >  th)),n);
end
[Nh,Vmh]=hist(Vm,0:4:300);
sm=spec1d(Vmh,Nh,sqrt(Nh));
[stmp,fm]=fits(sm,'gauss',[150 200 50 0],[1 1 1 0]);
[Nh,Vph]=hist(Vp,0:4:300);
sp=spec1d(Vph,Nh,sqrt(Nh));
[stmp,fp]=fits(sp,'gauss',[150 200 50 0],[1 1 1 0]);
clf
hist(Vm,[0:4:300])
axis([0 300 0 200])
nicefig
set(gca,'position',[.2 .15 .7 .4])
line(Vmh,gauss(Vmh,fm.pvals),'linewidth',1,'color','r')
text(14,175,sprintf('V_-= %3.0f \\pm %3.0f',[fm.pvals(2:3)]),'fontname','times','fontsize',12)
text(14,150,sprintf('Total %4.4i',length(Vm)),'fontname','times','fontsize',12)
text(14,125,sprintf('In Gaussian   %2.2i %%',round(sqrt(2*pi)*fm.pvals(3)*fm.pvals(1)/4/length(Vm)*100)),'fontname','times','fontsize',12)
text(14,100,sprintf('Below 0.8*V_-   %1.1i %%',round(sum(Vm<fm.pvals(2)*0.8)/length(Vm)*100)),'fontname','times','fontsize',12)
axes('position',[.2 .58 .7 .4])
hold on
hist(Vp,[0:4:300])
hold off
axis([0 300 0 200])
set(gca,'fontname','times','fontsize',12,'XTickLabel','')
box on
line(Vph,gauss(Vph,fp.pvals),'linewidth',1,'color','r')
text(14,175,sprintf('V_+= %3.0f \\pm %3.0f',[fp.pvals(2:3)]),'fontname','times','fontsize',12)
text(14,150,sprintf('Total %4.4i',length(Vp)),'fontname','times','fontsize',12)
text(14,125,sprintf('In Gaussian   %2.2i %%',round(sqrt(2*pi)*fp.pvals(3)*fp.pvals(1)/4/length(Vm)*100)),'fontname','times','fontsize',12)
text(14,100,sprintf('Below 0.8*V_+   %1.1i %%',round(sum(Vp<fp.pvals(2)*0.8)/length(Vp)*100)),'fontname','times','fontsize',12)
title('LSMO 2003/3, T=220 K, Vz=0.8 V, Iz=0.2 nA, (100 Å)^2')
print -depsc lsmo03feb_220K_gaphist

