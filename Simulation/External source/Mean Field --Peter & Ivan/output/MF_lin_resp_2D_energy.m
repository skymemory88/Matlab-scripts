function MF_lin_resp_2D_energy
close all
cd('D:\Projects\LiHoF4_Network_Analyzer\plots\new_analysis\MF_lin_resp_2D_energy')
load ('susc_calc300mK.mat','imchi','rechi1')

fields= [0:0.01:9];
freqs = [0:0.01:6];
% imchi  (zzz,k);

%colorplot dB
% zq = F(xq,yq);
% zq1 = medfilt1(zq,10,size(zq,1),2);
fig1 = figure (321);

% hp = pcolor(fields,freqs,-real(log(imchi)));
% set(hp, 'edgeColor','none')

[C hp]=contourf(fields,freqs,real(log(imchi)),300);
set(hp,'lineStyle','none')

set(gca,'fontsize',15)
set(fig1,'position',[10 10 600 400])
colorbar('location','NorthOutside')
xlabel('Field (T)','FontSize',15)
ylabel('Frequency (GHz)','FontSize',15)
colormap jet
caxis([-15 0])
% xlim([0 7])
% ylim([3.4 3.47])
grid off
box on
tit1=['MF_2D_energy'];    hold on     
hpcfr2 = plot([0 9],[5.604 5.604],'Color',[0.5 0 1],'linestyle','--');
hpcfr2 = plot([0 9],[4.449 4.449],'Color',[0 1 0.5],'linestyle','--');
hpcfr2 = plot([0 9],[3.924 3.924],'Color',[1 0 0.5],'linestyle','--');
hpcfr2 = plot([0 9],[3.436 3.436],'Color',[0 0.5 1],'linestyle','--');
hpcfr2 = plot([0 9],[1.682 1.682],'Color',[1 0.5 0],'linestyle','--');
% title(tit1,'FontSize',15)
% hold on
% save([tit1,'.mat'],'xq','yq','zq')
% h_bar = colorbar('East');
% initpos = get(h_bar,'Position');
% set(h_bar, 'Position',[initpos(1) initpos(2) initpos(3)*0.5 initpos(4)*0.5])

saveplots(fig1,tit1)
end

function saveplots(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
% dir=['D:\Measurements\MATLAB\data\',num2str(freq_cut,'%3.3f'),'GHz ',num2str(T,'%3.3f'),'K'];
% mkdir(dir)
cd('D:\Projects\LiHoF4_Network_Analyzer\plots\new_analysis\MF_lin_resp_2D_energy')

saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg'});

disp(['Figure ' figname ' saved to '])
disp(cd)
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption')
end
