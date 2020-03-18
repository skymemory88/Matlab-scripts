function figure2a
% close all
% temp=[0.100, 0.300, 0.500, 0.800, 1.3];
temp = [0.200];
color = ["black","red","blue","magenta","green","yellow","cyan"];
theta = [90]; % Angle between the transverse field and a-axis
marker = [":","-.","--","-"];
% theta = 0;
% color = {[255 127 0], [255 0 127], [127 0 255], [0 127 255], [0 255 127]} ;
lg = strings(1,2); % Create an empty array for legends
figs = gobjects(7,numel(temp));

    for iter = 1:numel(temp)
        for iter2 = 1:numel(theta)
%             cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output\A=1.0_angles_Yikai');
            cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output\');
%             cd('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Simulation/External source/Mean Field --Peter & Ivan/output/A=1.0_angles_Peter')
%             lname=[num2str(temp(iter),'%3.3f'),num2str(theta(iter2),'_%u'),'.mat']; 
            lname=[num2str(temp(iter),'%3.3f'),'.mat']; 
            
            load(lname,'-mat','eee','fff','h_mf2');
            fields = vecnorm(fff);
            hmf=h_mf2(:,3);
            E2f = 241.8; % Convert Energy to frequency
            E(:,:)=squeeze(eee)*E2f;
            %E(:,:) = eee(:,1,:)*E2F;
            Ediff = double.empty(0,size(E,1));

            % Calcuate average and deviation of the energy difference along the
    %         E8=E(:,1:8);
    %         for j=1:901
    %         aver(j)=sum(E8(j,:))/8;
    %         end
    %         for j=1:8
    %         E8(:,j)=E8(:,j)-aver(:);
    %         end
%%
    %         % Plot the lowest eight energy levels
    %         figure
    %         hold on
    %         hpcfr1 = plot(fields,E(:,1:8),'Color',[35 107 142]/255,'linewidth',2);
    %         set(gca,'fontsize',15); 
    %         % text(8,(j-1)*D+0.02,'0.1');
    % %         set(fig1,'position',[100 100 600 300])
    %         xlabel('Field (T)','FontSize',15)
    %         ylabel('Energy (GHz)','FontSize',15)
    %         grid off
    %         box on; 
    %         xlim([0 9]);
    % %         ylim([2 5]);
    %         tit1='Energy levels';         
    %         title(tit1,'FontSize',15)
    %         legend(num2str(temp(iter),'T = %.2f K,  A = A_{th}'))
%%
            % Calculate the energy difference bewteen neighbour levels
            for i=1:7
                Ediff(i,:)=E(:,i+1)-E(:,i);
            end
    %         % Plot the lowest energy difference between the 8 levels
%             figure
            % frequency = { '1.682 GHz', '3.436 GHz', '3.924 GHz', '4.449 GHz', '5.604 GHz' } ;
%             figure
            figs(:, iter2) = plot(fields, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter2),'LineWidth',2);
            hold on
            lg(iter, iter2) = [num2str(temp(iter),"T = %.3f K, A = A_{th}"),num2str(theta(iter2),", theta = %u")];


%             hpcfr1 = plot(fields, Ediff,'Color','blue','linewidth',2);
    %         hpcfr1 = plot(fields, Ediff(2:end,:),'Color','black','linewidth',2);   
    %         hpcfr1 = plot(fields, Ediff,'Color',[35 107 142]/255,'linewidth',2,'LineStyle','-');
%             legend(num2str(temp(iter),'T = %.2f K, A = A_{th}'),num2str(temp(iter),'T = %.2f K, A = A_{th}'))

    %         hpcfr2 = plot([0 9],[5.604 5.604],'Color',[0.5 0 1],'linestyle','--');
    %         hpcfr2 = plot([0 9],[4.449 4.449],'Color',[0 1 0.5],'linestyle','--');
    %         hpcfr2 = plot([0 9],[3.924 3.924],'Color',[1 0 0.5],'linestyle','--');
%             Qhpcfr2 = plot([0 17],[2.37 2.37],'Color',[1 0.5 1],'linestyle','--');
    %         hpcfr2 = plot([0 9],[1.682 1.682],'Color',[1 0.5 0],'linestyle','--');

%             % Load experimental data of line crossings and superimpose the plot onto the calculated data
%             load('f0_BC.mat');
%             plot(B(:,1),B(:,2),'ro','MarkerFaceColor','red');
            
%             yyaxis(gca,'left');
            set(gca,'fontsize',15,'Xtick',0:1:17); 
%             text(8,(j-1)*D+0.02,'0.1');
    %         set(fig4,'position',[800 100 600 300]);
            xlabel('Field (T)','FontSize',15);
            ylabel('Energy difference (GHz)','FontSize',15);
            grid off
            box on;  
            xlim([min(fields) max(fields)]);
            ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
            tit4='Difference between energy levels';         
            title(tit4,'FontSize',15)
            hold on
%             % Calculate the frequency range that the transition energy covers at each field
%             Erange = double.empty(0,length(Ediff));
%             for i=1:size(Ediff,2)
%                 Erange(i) = (max(Ediff(:,i))-min(Ediff(:,i)));
%             end
%             % Plot the frequency range that the energy differences cover at each field
%     %         figure
%     %         hold on
%             yyaxis right;
%             plot(fields, Erange,'k-');
%             plot(fields(Erange == min(Erange)),min(Erange),'o','MarkerSize',4,'MarkerFaceColor','none');
%             set(gca,'fontsize',15); 
%             xlabel('Field (T)','FontSize',15)
%             ylabel('Frequency range (GHz)','FontSize',15)
%             grid off
%             box on;  
%             xlim([0 17]);
%             tit4='Frequency range the energy levels cover';         
%             title(tit4,'FontSize',15)
        end
    end
    lgd = legend(figs(1,:),lg);
    lgd.FontSize = 12;
end

function saveplots(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
cd('D:\Projects\LiHoF4_Network_Analyzer\plots\new_analysis\figure2a')
saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg'});
disp(['Figure ' figname ' saved to '])
disp(cd)
% cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption')
end
