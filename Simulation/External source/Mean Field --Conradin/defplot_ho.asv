close all
clear all

load werte_mithyp2

X=(temperatur'*ones(1,size(punkte,2)));
Y=punkte(:,:,1);
Z=Swerte(:,:,3);


set(gcf,'papertype','a4','paperposition',[0.25 0.25 2*8.5 6.5],'paperunits','centimeter')
axes('position',[.22/2 .2 .76/2 .76])
set(gca,'fontname','helvetica','fontsize',12);
contour(X',Y',Z');
axis([0 2.2 0 5.2])
shading interp
ylabel('H_{x}    [T]','interp','tex')
xlabel('T    [K]','interp','tex')
cbar_axes = colorbar
set(cbar_axes, 'OuterPosition', [0.92 0.2 0.07 0.7])
set(cbar_axes,'fontname','helvetica','fontsize',12);
ctle=get(cbar_axes,'title')
set(cbar_axes,'fontname','helvetica','fontsize',12);
set(ctle,'String','<S^z>','interp', 'tex')
set(ctle,'fontname','helvetica','fontsize',12);
hold on

%Vergleich mit Daten von Bitko:

%A = imread('bitko','png');
%image(A*20)
%xy=ginput
%xnull=xy(1,1);
%ynull=xy(1,2);
%skalax=(xy(3,1)-xnull)/1.6; %in Kelvin
%skalay=(xy(2,2)-ynull)/6; %in Tesla
%for n=4:size(xy,1)
%    hbitko(n-3)=(xy(n,2)-ynull)/skalay;
%    tbitko(n-3)=(xy(n,1)-xnull)/skalax;
%end
%save( 'bitko', 'hbitko', 'tbitko')
load( 'bitko', 'hbitko', 'tbitko')
hbitko(1:end-2)=hbitko(1:end-2)-0.24%Demagnetisierung
plot(tbitko,hbitko,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)


load werte_ohnehyp2

X=(temperatur'*ones(1,size(punkte,2)));
Y=punkte(:,:,1);
Z=Swerte(:,:,3);

axes('position',[.22/2 .2 .76/2 .76]+[0.4 0 0 0])
set(gca,'fontname','helvetica','fontsize',12);
contour(X',Y',Z');
axis([0 2.2 0 5.2])
shading interp
%ylabel('H_{x}    [T]','interp','tex')
set(gca, 'YTick', [])
xlabel('T    [K]','interp','tex')
hold on
plot(tbitko,hbitko,'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)

print -depsc2 mfplot.eps


%%%#############################################################
%%Standart plot

set(gcf,'papertype','a4','paperposition',[0.25 0.25 12 9],'paperunits','centimeter')
h=plot(spa);
set(gca,'position',[.15 .15 .8 .8])
axis([0.26 0.47 0 1.6])
set(gca,'fontname','helvetica','fontsize',12);
grid off
xlabel('T [K]')
ylabel('c/R')
legend(h,num2str(10*feld','%1.2f kOe'))
set(h,'Markersize', 5)

