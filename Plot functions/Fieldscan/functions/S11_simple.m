%% One curve resonant frequency simulation
% % delete(splot)
% % delete(lplot)
% hold on
% freq_l = 3.5;
% freq_h = 3.56;
% gc = 0.05;
% f_cav = 3.54;
% wmm = min(Ediff(:,:)); % choose the groud state as the lowest in energy
% wpp = max(Ediff(:,:)); % choose the highest occupied level
% wp = 0.5.*(wpp + f_cav)+0.5.*sqrt((wpp-f_cav).^2 + gc^2);
% wm = 0.5.*(wmm + f_cav)-0.5.*sqrt((wmm-f_cav).^2 + gc^2);
% lplot = plot(fields,Ediff(:,:),':k','linewidth',2);
% hold on
% splot = plot(fields,wp,'-r',fields,wm,'-r','linewidth',2);
% xlim([min(fields) max(fields)])
% xlabel('Magnetic field (T)')
% ylabel('Energy (GHz)')
%% Split curve resonant frequency simulation 
% % In the case of strong coupling, the coupling strength change before and
% % afterthe quantum transition point (QTP), therefore split the curve into
% % two partes before and after the QTP.
% % 
% % delete(splot)
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
% splot2 = plot(B_l,wm,'-r',B_h,wm2,'-b','Linewidth',2);
%% S11 simulation
% freq_l = 3.4;
% freq_h = 3.58;
% g_const = 0.05;
% gc = 1.5*g_const;
% f_cav = 3.53;
% [B,w] = meshgrid(linspace(min(fields),max(fields),length(Ediff(1,:))),linspace(freq_l,freq_h,300));
% f2E = 1/241.8; % Convert from GHz to meV
% kB = 0.08617; % [meV/K]
% Temp = 0.150; % [K]
% kc = 0.1*g_const;
% gamma = 0.4*g_const;
% sum = 0.0;
% for ii = 1:7
%     gcc = gc.*sqrt(exp(-Ediff(ii,:).*f2E./(kB*Temp)));
%     sum = sum + gcc.^2./(1i.*(w-Ediff(ii,:))-gamma); % Sum over all dispersion curves
% end
% S11 = abs(1+kc./(1i.*(w-f_cav)-kc+sum));
% figure
% % surf(B, w, log(S11),'edgecolor','none');
% map = pcolor(B, w, log(S11));
% map.EdgeColor = 'none';
% map.FaceColor = 'interp';
% xlim([min(fields) 9]);
% % xlim([min(fields) max(fields)]);
% ylim([freq_l freq_h]);
% caxis([-8 0])
% xlabel('Magnetic field (T)');
% ylabel('Frequency (GHz)');
% title('Simulated S11 response');
% colorbar
% % view(0,90);
% clearvars -except Ediff fields