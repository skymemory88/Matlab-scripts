[B,w] = meshgrid(linspace(min(H0),max(H0),500),linspace(min(f0),max(f0),500));
kc = 0.1*fitPara.g;
ke = kc;
gc = fitPara.g;
wc = fitPara.wc;
Br = fitPara.x0;
gamma = fitPara.gamma;
hbar = 1.0;
Delt = -7/2.*(B(1,:)-Br)./hbar;
wp = wc + 0.5.*Delt + 0.5.*sqrt(Delt.^2 + 4*gc^2);
wm = wc + 0.5.*Delt - 0.5.*sqrt(Delt.^2 + 4*gc^2);
hold on
p1 = plot(B(1,:),wp,'-r',B(1,:),wm,'-r','Linewidth',1);

figure
S11 = abs(1+ke./(1i.*(w-wc)-kc+gc^2./(1i.*(-Delt)-gamma)));
surf(B, w, log(S11),'edgecolor','none');
xlim([min(H0) max(H0)]);
ylim([min(f0) max(f0)]);
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title('Simulated S11 response');
colorbar
view(0,90);