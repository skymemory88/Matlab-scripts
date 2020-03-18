[ivst1,pvst1,fvst1,evst1]=collate3(126,146);
[ivst2,pvst2,fvst2,evst2]=collate3(156,198);

ivst=[ivst1;ivst2];
pvst=[pvst1;pvst2];
fvst=[fvst1;fvst2];
evst=[evst1;evst2];

figure
subplot(2,2,1)
errorbar(fvst(:,1),fvst(:,2),evst(:,2),'o');
set(gca,'Yscale','Linear');
xlabel('Temperature (K)')
ylabel('Intensity (a.u)')
subplot(2,2,2)
errorbar(fvst(:,1),2*fvst(:,4),2*evst(:,4),'o');
xlabel('Temperature (K)');
ylabel('FWHM');
set(gca,'Yscale','Linear');
axis([5 8 0 0.2])
subplot(2,2,3)
errorbar(fvst(:,1),fvst(:,3),evst(:,3),'o');
xlabel('Temperature (K)');
ylabel('Centre');
set(gca,'Yscale','Linear');subplot(2,2,2)
subplot(2,2,4);
errorbar(fvst(:,1),fvst(:,5),evst(:,5),'o');
xlabel('Temperature (K)');
ylabel('Background');
set(gca,'Yscale','Linear');