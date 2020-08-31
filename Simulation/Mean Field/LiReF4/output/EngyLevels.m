function EngyLevels
% close all
hold on
clearvars
% temp=[0.100, 0.300, 0.500, 0.800, 1.3];
temp = 0.240;

Options.Elevel = false;
Options.Ediff = true;
Options.Espan = false;
Options.savedata = false;

color = ["black","red","blue","magenta","green","yellow","cyan"];
theta = [0.0]; % Angle (in degrees) deviated from the transverse field direction
phi = [0.0]; % Angle (in degrees) in a-b plane
marker = [":","-.","--","-"];
% theta = 0;
% color = {[255 127 0], [255 0 127], [127 0 255], [0 127 255], [0 255 127]} ;
lg = strings(1,numel(theta),numel(phi)); % Create an empty array for legends
figs = gobjects(7,numel(theta),numel(phi));
for iter = 1:numel(temp)
    for iter2 = 1:numel(theta)
        for iter3 = 1:numel(phi)
            location = 'G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output\without Hz_I';
            %         lname=[num2str(temp(iter),'%3.3f'),num2str(theta(iter2),'_%u'),'.mat'];
            lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.1fDeg_%3$.1fDeg',temp(iter),theta(iter2),phi(iter3)),'.mat'];
            file = fullfile(location,lname);
            %         load(lname,'-mat','eee','fff');
            load(file,'-mat','eee','fff');
            
            fields = vecnorm(fff);
            E2f = 241.8; % Convert Energy to frequency
            E(:,:)=squeeze(eee)*E2f;
            %E(:,:) = eee(:,1,:)*E2F;
            Ediff = double.empty(0,size(E,1));
            %% Deviation of the energy from the mean value
            %             aver = double.empty(size(E,1),0);
            %             E8=E(:,1:8);
            %             for j=1:size(E,1)
            %                 aver(j)=mean(E8(j,:));
            %             end
            %             for j=1:8
            %                 E8(:,j)=E8(:,j)-aver(:);
            %             end
            %             clearvars aver
            %
            %             figure
            %             hold on
            %             plot(fields,E8(:,1:8),'Color',[35 107 142]/255,'linewidth',2);
            %             xlabel('Field (T)','FontSize',15)
            %             ylabel('Energy (GHz)','FontSize',15)
            %             xlim([0 9]);
            %             tit0='Energy spread from the mean value';
            %             title(tit0,'FontSize',15)
            %             legend(num2str(temp(iter)*1000,'T = %u mK,  A = A_{th}'))
            %% Plot the lowest eight energy levels
            if Options.Elevel == true
                figure
                hold on
                plot(fields,E(:,1:8),'Color',[35 107 142]/255,'linewidth',2);
                set(gca,'fontsize',15);
                % text(8,(j-1)*D+0.02,'0.1');
                %         set(fig1,'position',[100 100 600 300])
                xlabel('Field (T)','FontSize',15)
                ylabel('Energy (GHz)','FontSize',15)
                grid off
                box on;
                xlim([0 9]);
                %         ylim([2 5]);
                tit1='Energy levels';
                title(tit1,'FontSize',15)
                legend(num2str(temp(iter)*1000,'T = %u mK,  A = A_{th}'))
            end
            %% Energy difference bewteen neighbour levels
            if Options.Ediff == true
                for i=1:7
                    Ediff(i,:)=E(:,i+1)-E(:,i);
                end
                %         % Plot the lowest energy difference between the 8 levels
                %             figure
                % frequency = { '1.682 GHz', '3.436 GHz', '3.924 GHz', '4.449 GHz', '5.604 GHz' } ;
                figs(:, iter2, iter3) = plot(fields, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter3),'LineWidth',2);
                hold on
                lg(iter, iter2, iter3) = [num2str(temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %.1f"),char(176),num2str(phi(iter3),", phi = %.1f"),char(176)];
                
                %             hpcfr1 = plot(fields, Ediff,'Color','blue','linewidth',2);
                %             hpcfr1 = plot(fields, Ediff(2:end,:),'Color','black','linewidth',2);
                %             hpcfr1 = plot(fields, Ediff,'Color',[35 107 142]/255,'linewidth',2,'LineStyle','-');
                %             legend(num2str(temp(iter),'T = %.2f K, A = A_{th}'),num2str(temp(iter),'T = %.2f K, A = A_{th}'))
                %
                %             % Plot the resonant frequency of a bare cavity
                %             f_cav = 3.54; % Set resonant frequency of the bare cavity
                %             plot([0 max(fields)],[f_cav f_cav],'-r','LineWidth',1.5);
                %
                %             % Load experimental data of line crossings and superimpose the plot onto the calculated data
                %             load('f0_BC.mat');
                %             plot(B(:,1),B(:,2),'ro','MarkerFaceColor','red');
                %
                %             yyaxis(gca,'left');
                set(gca,'fontsize',15,'Xtick',0:1:max(fields));
                %             text(8,(j-1)*D+0.02,'0.1');
                %             set(fig4,'position',[800 100 600 300]);
                xlabel('Field (T)','FontSize',15);
                ylabel('Energy difference (GHz)','FontSize',15);
                grid off
                box on;
                xlim([min(fields) max(fields)]);
                ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
                tit4='Difference between energy levels';
                title(tit4,'FontSize',15)
            end
        end
        %% Frequency span of the transitions at each field
        if Options.Espan == true
            hold on
            eRange = double.empty(0,length(Ediff));
            for i=1:size(Ediff,2)
                eRange(i) = (max(Ediff(:,i))-min(Ediff(:,i)));
            end
            yyaxis right;
            plot(fields, eRange,'k-');
            plot(fields(eRange == min(eRange)),min(eRange),'o','MarkerSize',4,'MarkerFaceColor','none');
            set(gca,'fontsize',15);
            xlabel('Field (T)','FontSize',15)
            ylabel('Frequency range (GHz)','FontSize',15)
            grid off
            box on;
            xlim([min(fields) max(fields)]);
            ylim([0.9*min(eRange) 1.2*max(eRange)]);
            tit4='Frequency range the energy levels cover';
            title(tit4,'FontSize',15)
        end
    end
end
lgd = legend(figs(1,:,:),lg);
lgd.FontSize = 12;
if Options.savedata == true
    filename = strcat('sim_',num2str(temp*1000,'%u'),'mK_trans.mat','-v7.3');
    save(filename,'fields','Ediff');
end
end
