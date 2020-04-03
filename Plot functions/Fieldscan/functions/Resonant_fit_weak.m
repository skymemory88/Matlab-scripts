hold on
B = linspace(min(H0),max(H0),200);
kc = 0.1*fitPara.g;
ke = kc;
gc = fitPara.g;
wc = fitPara.wc;
Br = fitPara.x0;
gamma = fitPara.gamma;
hbar = 1.0;
Delt = 7/2.*(B-Br)./hbar;
w = wc + gc^2*Delt./(Delt.^2+gamma^2);
plot(B,w,'-r','LineWidth',2);
xlim([min(H0) max(H0)]);
ylim([min(min(f0),min(w)) 4.085]);
xlabel('Magnetic field (T)');
ylabel('Frequency (GHz)');
title('Simulated S11 response');