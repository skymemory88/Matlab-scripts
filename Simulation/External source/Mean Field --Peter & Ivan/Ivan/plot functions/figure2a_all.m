function figure2a_all
close all
temp=[0.001, 0.05:0.05:2 , 2.5, 3, 4];

for j=1: length(temp)
T=temp(j)
cd('W:\MF_calc\dataA097\data')
lname=[num2str(temp(j),'%3.3f'),'.mat']; 
load(lname)

fields=fff(1,:);
hmf=h_mf2(:,3);
ghztomeV = 1/241.8;
E(:,:)=eee(:,1,:)/ghztomeV;

for i=1:7
    Ediff(i,:)=E(:,i+1)-E(:,i);
end

E8=E(:,1:8);
% for j=1:901
% aver(j)=sum(E8(j,:))/8;
% end
% for j=1:8
% E8(:,j)=E8(:,j)-aver(:);
% end

cd('W:\MF_calc\dataA097\energies')
save (lname,'T','fields','Ediff','E8','E')
end

% fig1 = figure (3239878); hold on;
% hpcfr1 = plot(fields,E8,'Color',[35 107 142]/255,'linewidth',2);
% set(gca,'fontsize',15); 
% % text(8,(j-1)*D+0.02,'0.1');
% set(fig1,'position',[100 100 600 300])
% xlabel('Field (T)','FontSize',15)
% ylabel('Energy (GHz)','FontSize',15)
% grid off
% box on; 
% xlim([0 9]);
% % ylim([0 5.5]);
% tit1=['Energy levels'];         
% title(tit1,'FontSize',15)
% % saveplots(fig4,tit4,frequency{j},'0.1')
% % close
% legend('T=0.15K, A=0.97')
% saveplots(fig1,tit1)
% 
% fig4 = figure (32378); hold on;
% 
% % Ccc = {'orange','red','purple','blue','green'} ;
% % C = {[255 127 0], [255 0 127], [127 0 255], [0 127 255], [0 255 127]} ;
% % frequency = { '1.682 GHz', '3.436 GHz', '3.924 GHz', '4.449 GHz', '5.604 GHz' } ;
% 
% hpcfr1 = plot(fields,Ediff  ,'Color',[35 107 142]/255,'linewidth',2);
% hpcfr2 = plot([0 9],[5.604 5.604],'Color',[0.5 0 1],'linestyle','--');
% hpcfr2 = plot([0 9],[4.449 4.449],'Color',[0 1 0.5],'linestyle','--');
% hpcfr2 = plot([0 9],[3.924 3.924],'Color',[1 0 0.5],'linestyle','--');
% hpcfr2 = plot([0 9],[3.436 3.436],'Color',[0 0.5 1],'linestyle','--');
% hpcfr2 = plot([0 9],[1.682 1.682],'Color',[1 0.5 0],'linestyle','--');
% set(gca,'fontsize',15); 
% % text(8,(j-1)*D+0.02,'0.1');
% set(fig4,'position',[800 100 600 300])
% xlabel('Field (T)','FontSize',15)
% ylabel('Energy difference (GHz)','FontSize',15)
% grid off
% box on; 
% xlim([0 9]);
% ylim([0 7]);
% tit4=['Difference between energy levels'];         
% title(tit4,'FontSize',15)
% % saveplots(fig4,tit4,frequency{j},'0.1')
% % close
% legend('T=0.300K, A=0.97')
% saveplots(fig4,tit4)
end


function saveplots(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
% cd('W:\MATLAB\data\3.924GHz 0.010K')
cd('D:\Projects\LiHoF4_Network_Analyzer\plots\new_analysis\figure2a300')
saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg'});
disp(['Figure ' figname ' saved to '])
disp(cd)
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption')
end
