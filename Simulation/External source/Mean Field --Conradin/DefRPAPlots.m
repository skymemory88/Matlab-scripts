%ohne Hyperfeinww q=(200) T=0.05 epsilon=0.01
clear all
close all

load('disph_ohnehyp','hvek','ow','Skw');

set(gcf,'papertype','a4','paperposition',[0.25 0.25 2*8.5 6.5],'paperunits','centimeter')
axes('position',[.22/2 .2 .76/2 .76])
set(gca,'fontname','helvetica','fontsize',12);


pcolor(hvek,ow,log(Skw)); 
shading interp
set(gca,'YLim',[0.01,0.49]);
xlabel('H_x / Tesla','interp','tex')
ylabel('\omega / meV','interp','tex')
hold on

%mit Hyperfeinww q=(200) T=0.05 epsilon=0.01
load('disph_mithyp','hvek','ow','Skw');

axes('position',[.22/2 .2 .76/2 .76]+[0.5 0 0 0])
set(gca,'fontname','helvetica','fontsize',12);

pcolor(hvek,ow,log(Skw)); 
shading interp
set(gca,'YLim',[0.01,0.49]);
xlabel('H_x / Tesla','interp','tex')
ylabel('\omega / meV','interp','tex')

print -depsc2 rpaplotho.eps
