B = linspace(0,9,100);
Br = 3.72;
m = 7/2;
hbar = 1;
Delt = -m*(B-Br)/hbar;
wc = 3.674;
gc = 0.055;
gamma = 0.1;
w = wc - gc^2.*Delt./(Delt.^2+gamma^2);
hp1 = plot(B,w,'-k','LineWidth',3);
hold on
wp = wc + Delt./2 + sqrt(Delt.^2+4*gc^2)/2;
wm = wc + Delt./2 - sqrt(Delt.^2+4*gc^2)/2;
hp2 = plot(B,wm,'-r',B,wp,'-r','LineWidth',3);