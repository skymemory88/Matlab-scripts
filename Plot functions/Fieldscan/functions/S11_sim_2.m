delete(splot)
hold on
% gc = 0.3;
gp = 0.1;
gm = 0.3;
QPT = 4.9;
f_cav = 3.53;
B_l = fields(fields <= QPT);
B_h = fields(fields > QPT);
% wmm = min(Ediff(:,:));
% wpp = Ediff(:,:);
% wpp = Ediff(6,1:numel(B_l));
% wpp2 = Ediff(6,numel(B_l)+1:end);
wmm = Ediff(1,1:numel(B_l));
wmm2 = Ediff(6,numel(B_l)+1:end);
% wp = 0.5.*(wpp + f_cav)+0.5.*sqrt((wpp-f_cav).^2 + gm^2);
% wp2 = 0.5.*(wpp2 + f_cav)+0.5.*sqrt((wpp2-f_cav).^2 + gp^2);
wm = 0.5.*(wmm + f_cav)-0.5.*sqrt((wmm-f_cav).^2 + gm^2);
wm2 = 0.5.*(wmm2 + f_cav)-0.5.*sqrt((wmm2-f_cav).^2 + gp^2);
% level = plot(fields,Ediff,'-k');
% splot = plot(fields,wp,'-r',fields,wm,'-r','linewidth',2);
% splot = plot(B_l,wm,'-r',B_h,wm2,'-b',B_l,wp,'-k',B_h,wp2,'-k','Linewidth',2);
splot = plot(B_l,wm,'-r',B_h,wm2,'-r','Linewidth',2);