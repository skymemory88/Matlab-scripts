function EngyLevels
clearvars

hbar = 1.055E-34; % Reduced Planck constant [J.s]
meV2J = 1.602217e-22; % [J/meV]
f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]

% temp =  [1.774 1.776 1.778 1.780 1.781 1.782 1.783 1.784 1.785 1.787];
temp = 0.16;
theta = 0.0; % Angle (in degrees) between magnetic field and c-axis in a-c plane
phi = 0.0; % Angle (in degrees) in a-b plane
N_level = 7; % Number of transition levels to plot

Options.ion = 'Ho'; % support 'Er' and 'Ho'
Options.hyp = 1.00; % Hyperfine isotope proportion
Options.nZee = false; % nuclear Zeeman interaction
Options.RPA = false; % not enabled yet
Options.Js = false; % spin expectation values
Options.eUnit = 'GHz';
    if strcmp(Options.eUnit,'meV')
        E2f = 1; % use default meV unit
        ylab = 'Energy (meV)';
        ylab1 = 'Energy gap (meV)';
        ylab2 = 'Freqeuency range (meV)';
    elseif strcmp(Options.eUnit,'GHz')
        E2f = 1/f2E; % meV to GHz conversion
        ylab = 'Energy (GHz)';
        ylab1 = 'Energy gap (GHz)';
        ylab2 = 'Freqeuency range (GHz)';
    end
Options.Elevel = false; % eigen-energies
Options.Ediff = true; % excitation spectrum
    Options.deltaI2 = 0; % Transitions between the next nearest neightbouring levels
Options.Espan = 0;
Options.savedata = 0;
Options.savegif = 0;
if Options.nZee == true
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',Options.ion,...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I'];
else
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',Options.ion,...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\without Hz_I'];
end

marker = [":","-.","--","-"];
color = ["black","red","blue","magenta","green","yellow","cyan"];
% color = {[255 127 0], [255 0 127], [127 0 255], [0 127 255], [0 255 127]};
lg = strings(1,numel(theta),numel(phi)); % Create an empty array for legends
figs_E = gobjects(N_level+1,numel(theta),numel(phi));
figs_dE = gobjects(N_level,numel(theta),numel(phi));
if Options.deltaI2; figs_dE2 = gobjects(N_level-1,numel(theta),numel(phi)); end
if Options.savegif == true
    im_t = numel(temp)*numel(theta)*numel(phi);
    im_idx = 1;
    if Options.Elevel == true
        Elevel_frame = []*im_t;
    end
    if Options.Ediff == true
        Ediff_frame = []*im_t;
    end
end
for iter = 1:numel(temp)
    for iter2 = 1:numel(theta)
        for iter3 = 1:numel(phi)
            if Options.RPA == true
                lname=['Hscan_Li',Options.ion,'F4_',...
                    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f',temp(iter),theta(iter2),phi(iter3),Options.hyp),'.mat'];
            else
                lname=['Hscan_Li',Options.ion,'F4_',...
                    sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f',temp(iter),theta(iter2),phi(iter3),Options.hyp),'.mat'];
            end
            file = fullfile(Options.location,lname);
            load(file,'-mat','vvv','ttt','eee','fff','ion');
            eee = eee- min(eee,[],2); % Normalize the energy amplitude to the lowest eigen-energy
            fields = vecnorm(fff);
            E(:,:)=squeeze(eee)*E2f;
            %E(:,:) = eee(:,1,:)*E2F;
            Ediff = double.empty(0,size(E,1));
            if Options.deltaI2; Ediff2 = double.empty(0,size(E,1)); end
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
%             ylabel(ylab,'FontSize',15)
%             xlim([0 9]);
%             tit0='Energy spread from the mean value';
%             title(tit0,'FontSize',15)
%             legend(num2str(temp(iter)*1000,'T = %u mK,  A = A_{th}'))
            %% Plot eigen-energies
            if Options.Elevel == true
%                 fig_E = figure;
%                 figs_E(:, iter2, iter3) = plot(fields,E(:,1:N_level+1),'Color',[35 107 142]/255,'linewidth',2);
                figs_E(:, iter2, iter3) = plot(fields,E(:,1:N_level+1),'Color','r','linewidth',2);
                hold on
                set(gca,'fontsize',15);
%                text(8,(j-1)*D+0.02,'0.1');
%                set(fig1,'position',[100 100 600 300])
                xlabel('Field (T)','FontSize',15)
                ylabel(ylab,'FontSize',15)
                grid off
                box on;
                tit1='Energy levels';
                title(tit1,'FontSize',15)
                legend(num2str(temp(iter)*1000,'T = %u mK,  A = A_{th}'))
                if Options.savegif == true
                    drawnow
                    frame = getframe(fig_E);
                    Elevel_frame{im_idx} = frame2im(frame);
                end
            end
            %% Plot excitation spectrum among neighbour levels
            if Options.Ediff == true
                for i=1:N_level-1 % Up until the second from top level
                    Ediff(i,:) = E(:,i+1)-E(:,i); % Transition between the nearest neighbouring levels
                    if Options.deltaI2
                        Ediff2(i,:) = E(:,i+2)-E(:,i); % Transition between the next nearest neighbouring levels
                    end
                end
                Ediff(N_level,:) = E(:,N_level+1)-E(:,N_level); % Transition between the top two levels
                
                % frequency = { '1.682 GHz', '3.436 GHz', '3.924 GHz', '4.449 GHz', '5.604 GHz' }; % frequency reference lines

                % Plot transitions between the nearest levels
%                 figure;
                hold on
                figs_dE(1, iter2, iter3) = plot(fields, Ediff(1,:), 'Marker', 'none', 'LineStyle', marker(1), 'Color', 'r','LineWidth',1.5);
                if N_level > 1
%                     figs_dE(2:end, iter2, iter3) = plot(fields, Ediff(2:end,:), 'Marker', 'none', 'LineStyle', marker(1), 'Color', color(iter3),'LineWidth',2);
                    figs_dE(2:end, iter2, iter3) = plot(fields, Ediff(2:end,:), 'Marker', 'none', 'LineStyle', marker(1), 'Color','k','LineWidth',1.5);
%                     figs_dE(2:end, iter2, iter3) = plot(fields, Ediff(2:end,:), 'Marker', 'none', 'LineStyle', marker(1),'LineWidth',1.5);
                end
                
%                 % Plot transitions between the next nearest levels
                if Options.deltaI2
                    figs_dE2(1, iter2, iter3) = plot(fields, Ediff2(1,:), 'Marker', 'none', 'LineStyle', marker(1), 'Color', 'r','LineWidth',2);
                    if N_level > 2
                        figs_dE2(2:end, iter2, iter3) = plot(fields, Ediff2(2:end,:), 'Marker', 'none', 'LineStyle', marker(1), 'Color', color(iter3),'LineWidth',2);
                    end

                    hpcfr1 = plot(fields, Ediff,'Color','blue','linewidth',2);
%                     hpcfr1 = plot(fields, Ediff(2:end,:),'Color','black','linewidth',2);
%                     hpcfr1 = plot(fields, Ediff,'Color',[35 107 142]/255,'linewidth',2,'LineStyle','-');
                    legend(num2str(temp(iter),'T = %.2f K, A = A_{th}'),num2str(temp(iter),'T = %.2f K, A = A_{th}'))

                    % Plot the resonant frequency of a bare cavity
                    f_cav = 3.54; % Set resonant frequency of the bare cavity
                    plot([0 max(fields)],[f_cav f_cav],'-r','LineWidth',1.5);

                    yyaxis(gca,'left');
                    xticks(linspace(0.1,max(fields),10))
%                     text(8,(j-1)*D+0.02,'0.1');
                    set(figs_dE2,'position',[800 100 600 300]);
                end
                set(gca,'fontsize',15,'Xtick',0:1:max(fields));
                set(gca,'XTickLabelRotation',0)
                xlabel('Magnetic Field (T)','FontSize',15);
                ylabel(ylab1,'FontSize',15);
                grid off
                box on;
                xlim([min(fields) max(fields)]);
                xticks('auto')
                ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
%                 ylim([1 5]);
%                 ylim([3.59 3.71]) % adjust y-axis range to fit experimental data
                tit4 = 'Energy gaps at ';
                title([tit4 sprintf('T = %.3f K',temp(iter))],'FontSize',15)
                if Options.savegif == true
                    drawnow
                    frame = getframe(fig_dE);
                    Ediff_frame{im_idx} = frame2im(frame);
                end
                lgd = legend(squeeze(figs_dE(1,:,:)),lg);
                lgd.FontSize = 12;
            end
            %% Add frequency span of the transitions at each field
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
                ylabel(ylab2,'FontSize',15)
                grid off
                box on;
                xlim([min(fields) max(fields)]);
                ylim([0.9*min(eRange) 1.2*max(eRange)]);
                tit4='Frequency range the energy levels cover';
                title(tit4,'FontSize',15)
            end
%             lg(iter, iter2, iter3) = [num2str(temp(iter)*1000,"T = %u mK, A = A_{th}"),', \theta', num2str(theta(iter2)," = %.1f"),...
%                                       char(176),', \phi', num2str(phi(iter3)," = %.1f"),char(176)];
            lg(iter, iter2, iter3) = ['\phi', num2str(phi(iter3)," = %.1f"), char(176)];
            if Options.savegif == true
                im_idx = im_idx + 1;
            end
            %% Plot expectation values of spin moment
            if Options.Js == true
                Jx = cell(1,4);
                Jy = cell(1,4);
                Jz = cell(1,4);
                Jnorm = cell(1,4);
                
                figure;
                hold on
                for ii = 1:4
                    Jx{ii} = double.empty(4,length(fff(1,:)),0);
                    Jy{ii} = double.empty(4,length(fff(1,:)),0);
                    Jz{ii} = double.empty(4,length(fff(1,:)),0);
                    for jj = 1:length(fff(1,:))
                        Jx{ii}(:,jj,1) = ion.Js_hyp(ii,1,jj,1);
                        Jy{ii}(:,jj,1) = ion.Js_hyp(ii,2,jj,1);
                        Jz{ii}(:,jj,1) = ion.Js_hyp(ii,3,jj,1);
                    end
                    Jnorm{ii} = vecnorm([Jx{ii};Jy{ii};Jz{ii}]);
                    plot(vecnorm(fff),Jnorm{ii},'o')
                end
                legend(strcat(Options.ion,'-1'),strcat(Options.ion,'-2'),strcat(Options.ion,'-3'),strcat(Options.ion,'-4'))
                xlabel('Magnetic field (T)')
                ylabel('<J>')
                
                Jmx = double.empty(4,length(fff(1,:)),0);
                Jmy = double.empty(4,length(fff(1,:)),0);
                Jmz = double.empty(4,length(fff(1,:)),0);
                for ii = 1:length(fff(1,:))
                    Jmx(:,ii,1) = ion.Jmom_hyp(1,ii,1);
                    Jmy(:,ii,1) = ion.Jmom_hyp(2,ii,1);
                    Jmz(:,ii,1) = ion.Jmom_hyp(3,ii,1);
                end
                figure
                hold on
                plot(vecnorm(fff),Jmx,'-s')
                plot(vecnorm(fff),Jmy,'-s')
                plot(vecnorm(fff),Jmz,'-s')
                legend('<Jx>','<Jx>','<Jz>')
                xlabel('Magnetic field (T)')
                ylabel('<J_i>')
            end
        end
    end
end

if Options.savedata == true
    filename = strcat('sim_',num2str(temp*1000,'%u'),'mK_trans.mat','-v7.3');
    save(filename,'fields','Ediff');
end
if Options.savegif == true
    if Options.Elevel == true
        filename = strcat('sim_',num2str(min(temp)*1000,'%u'),'-',num2str(max(temp)*1000,'%u'),'mK_Elevel.gif');
        fileobj = fullfile(Options.location,filename);
        for ii = 1:length(Elevel_frame)
            [img,cmp] = rgb2ind(Elevel_frame{ii},256);
            if ii == 1
                imwrite(img,cmp,fileobj,'gif','LoopCount',inf,'DelayTime',1);
            else
                imwrite(img,cmp,fileobj,'gif','WriteMode','append','DelayTime',1);
            end
        end
    end
    if Options.Ediff == true
        filename = strcat('sim_',num2str(min(temp)*1000,'%u'),'-',num2str(max(temp)*1000,'%u'),'mK_Ediff.gif');
        fileobj = fullfile(Options.location,filename);
        for ii = 1:length(Ediff_frame)
            [img,cmp] = rgb2ind(Ediff_frame{ii},256);
            if ii == 1
                imwrite(img,cmp,fileobj,'gif','LoopCount',inf,'DelayTime',1);
            else
                imwrite(img,cmp,fileobj,'gif','WriteMode','append','DelayTime',1);
            end
        end
    end
end
end
