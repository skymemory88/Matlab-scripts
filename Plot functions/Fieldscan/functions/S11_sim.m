% % delete(splot)
% % delete(lplot)
hold on
gc = 0.055;
f_cav = 3.54;
wmm = min(Ediff(:,:)); % choose the groud state as the lowest in energy
wpp = max(Ediff(:,:)); % choose the highest occupied level
wp = 0.5.*(wpp + f_cav)+0.5.*sqrt((wpp-f_cav).^2 + gc^2);
wm = 0.5.*(wmm + f_cav)-0.5.*sqrt((wmm-f_cav).^2 + gc^2);
lplot = plot(fields,Ediff(:,:),':k','linewidth',2);
hold on
splot = plot(fields,wp,'-r',fields,wm,'-r','linewidth',2);
xlim([min(fields) max(fields)])
xlabel('Magnetic field (T)')
ylabel('Energy (GHz)')
%% Two parts simulation 
% % In the case of strong coupling, the coupling strength change before and
% % afterthe quantum transition point (QTP), therefore split the curve into
% % two partes before and after the QTP.
% % 
% delete(splot)
% hold on
% f_cav = 3.57;
% gc_f = 0.4;
% gc_p = 0.25;
% QTP = 4.9;
% B_l = fields(fields <= QTP);
% B_h = fields(fields > QTP);
% % wpp = Ediff(1,1:numel(B_l));
% % wpp2 = Ediff(7,numel(B_l)+1:end);
% wmm = Ediff(1,1:numel(B_l));
% wmm2 = Ediff(7,numel(B_l)+1:end);
% % wp = 0.5.*(wpp + f_cav)+0.5.*sqrt((wpp-f_cav).^2 + gc_f^2);
% % wp2 = 0.5.*(wpp2 + f_cav)+0.5.*sqrt((wpp2-f_cav).^2 + gc_p^2);
% wm = 0.5.*(wmm + f_cav)-0.5.*sqrt((wmm-f_cav).^2 + gc_f^2);
% wm2 = 0.5.*(wmm2 + f_cav)-0.5.*sqrt((wmm2-f_cav).^2 + gc_p^2);
% % lplot = plot(fields,Ediff,'-k');
% % splot = plot(B_l,wm,'-r',B_h,wm2,'-b',B_l,wp,'-k',B_h,wp2,'-k','Linewidth',2);
% splot = plot(B_l,wm,'-r',B_h,wm2,'-b','Linewidth',2);