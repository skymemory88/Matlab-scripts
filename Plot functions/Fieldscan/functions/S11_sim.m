[B,w] = meshgrid(linspace(min(H0),max(H0),500),linspace(min(f0),max(f0),500));
kc = 0.005;
ke = 0.005;
gc = fitPara.g;
wc = fitPara.wc;
Br = fitPara.x0;
gamma = fitPara.gamma;
hbar = 1.0;
Delt = 7/2.*(B-Br)./hbar;
S11 = abs(1+ke./(1i.*(w-wc)-kc+gc^2./(1i.*Delt-gamma)));
surf(B, w, log(S11),'edgecolor','none');
xlim([min(H0) max(H0)]);
ylim([min(f0) max(f0)]);
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title('Simulated S11 response');
colorbar
view(0,90);