function [results, fields] = EngyLevels(ion, temp, theta, phi, N_level)
clearvars -except ion temp theta phi N_level

kB = 8.617333e-2; % [meV/K]
hbar = 1.055E-34; % Reduced Planck constant [J.s]
J2meV = 6.24151e+21; % [meV/J]
Gh2mV = hbar*2*pi*10^9*J2meV; % [meV/GHz]

% temp =  [1.774 1.776 1.778 1.780 1.781 1.782 1.783 1.784 1.785 1.787];
Options.ion = ion; % support 'Er' and 'Ho'
    elem_idx = find(strcmp(Options.ion,[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}]));
Options.hyp = 1.00; % Hyperfine isotope proportion
Options.nZee = false; % nuclear Zeeman interaction
Options.RPA = false; % not enabled yet
Options.Elevel = false; % eigen-energies
Options.Ediff = false; % excitation spectrum
    Options.deltaI2 = 0; % Transitions between the next nearest neightbouring levels
Options.Js = false; % spin expectation values
Options.Mx = true; % nearest off diagonal element of spin matrix
Options.Espan = false;
    f_cav = 3.645; % Resonant frequency of the cavity (GHz)
    filFctr = 0.0127*0.15; % filling factor
Options.eUnit = 'GHz';
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
if strcmp(Options.eUnit,'meV')
    E2f = 1; % use default meV unit
    ylab = 'Energy (meV)';
    ylab1 = 'Energy gap (meV)';
    ylab2 = 'Freqeuency range (meV)';
elseif strcmp(Options.eUnit,'GHz')
    E2f = 1/Gh2mV; % meV to GHz conversion
    ylab = 'Energy (GHz)';
    ylab1 = 'Energy gap (GHz)';
    ylab2 = 'Freqeuency range (GHz)';
end

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
results.En = {};
results.Ediff = {};
results.Js = {};
results.Jmx = {};
results.JImx = {};
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
            eigenE = squeeze(eee) * E2f;
            eigenW = vvv;
            fields = vecnorm(fff);
            Ediff = double.empty(0,size(eigenE,1));
            if Options.deltaI2; Ediff2 = double.empty(0,size(eigenE,1)); end
%%          Deviation of the energy from the mean value
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
                fig_E = figure;
%                 figs_E(:, iter2, iter3) = plot(fields,E(:,1:N_level+1),'Color',[35 107 142]/255,'linewidth',2);
                results.En{iter,iter2,iter3} = eigenE(:,1:N_level+1);
                figs_E(:, iter2, iter3) = plot(fields, eigenE(:,1:N_level+1),'Color','r','linewidth',2);
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
                    Ediff(i,:) = eigenE(:,i+1)-eigenE(:,i); % Transition between the nearest neighbouring levels
                    if Options.deltaI2
                        Ediff2(i,:) = eigenE(:,i+2)-eigenE(:,i); % Transition between the next nearest neighbouring levels
                    end
                end
                Ediff(N_level,:) = eigenE(:,N_level+1)-eigenE(:,N_level); % Transition between the top two levels
                results.Ediff{iter,iter2,iter3} = Ediff;
                
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
                results.Js{iter,iter2,iter3} = cat(3,Jmx,Jmy,Jmz);
                figure
                hold on
                plot(vecnorm(fff),Jmx,'-s')
                plot(vecnorm(fff),Jmy,'-s')
                plot(vecnorm(fff),Jmz,'-s')
                legend('<Jx>','<Jx>','<Jz>')
                xlabel('Magnetic field (T)')
                ylabel('<J_i>')
            end
%% Plot the matrix elements of the CMP coupling strength
            if Options.Mx == true
%                         Er      Ho      Yb      Tm      Gd      Y
                J_tab = [15/2;    8;      7/2;    6;      7/2;    1];
                L_tab = [6;       6;       3;     5;      0;      1];
                S_tab = [3/2;     2;      1/2;    1;      7/2;    1];
                I_tab = [3;      3.5;      0;     0;      1.50;   0];
                nLande = [-0.1611; 1.192; -0.2592; -0.462; -0.2265; 0]; % https://easyspin.org/documentation/isotopetable.html

                %Ions' lattice parameters
                lattice = [{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
                           {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
                           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
                           {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
                           {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
%                            {[5.132 0 0; 0 5.132 0; 0 0 10.59]}     %Gd
                           {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y
                
                lattice = lattice{elem_idx};
                ionJ = J_tab(elem_idx); % Electronic moment for Ho3+
                ionI = I_tab(elem_idx); % Nuclear moment for Ho3+

                muN = 5.05078e-27; % Nuclear magneton [J/T]
                muB = 9.274e-24; %[J/T]
                mu0 = 4*pi*1e-7; % [H/m]
                ELEf = gLande(L_tab(elem_idx),S_tab(elem_idx)) * muB; % Lande factor * Bohr magneton (J/T)
                NUCf = nLande(elem_idx) * muN; % (J/T)
                rho = 4e30 / det(lattice); % magnetic ion number density [m^-3]
                gw0 = sqrt(mu0*2*pi*f_cav*10^9*rho*filFctr/hbar/2) / (2*pi) * 1e-9; % susceptibility prefactor [T/J. Hz]*E-9

                %Initiate ionJ operators
                Jz=diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
                JhT.z=kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
                Jp=diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
                Jm=Jp'; % electronic spin ladder operator
                Jph=kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
                Jmh=kron(Jm,eye(2*ionI+1));
                JhT.x=(Jph+Jmh)/2;
                JhT.y=(Jph-Jmh)/2i;
                
                %Initiate I operators
                Iz=diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
                IhT.z=kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
                Ip=diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
                Im=Ip'; % Nuclear spin ladder operator
                Iph=kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
                Imh=kron(eye(2*ionJ+1),Im);
                IhT.x=(Iph+Imh)/2;
                IhT.y=(Iph-Imh)/2i;

                temperature = repmat(temp,length(fields),1);
                tz = double.empty(size(eigenW,2),size(eigenW,2),size(fields,2),0); % Expectation value of J pseudo-spin
                ttz = double.empty(size(eigenW,2),size(eigenW,2),size(fields,2),0); % Expectation value of J-I pseudo-spin
                Gc2 = double.empty(size(eigenW,2),size(eigenW,2),size(fields,2),0); % Coupling strength
                for kk = 1:size(fields,2) % calculate susceptibility for all fields
                    en = squeeze(eigenE(kk,:)); % Obtain the corresponding eigen energies [meV]
                    if temperature(kk) ~= 0
                        beta = 1/E2f / (kB*temperature(kk)); % [1/meV or 1/GHz]
                        Z = sum(exp(-beta*en));
                        zn = exp(-beta*en)/Z;
%                         Z = exp(-beta*en);
%                         for nn = 1:8 % Skew the population of the lowest 8 states
%                             Z(nn) = Z(nn) - (9-nn)/2e2 + nn/2e2;
%                         end
%                         Z = Z/sum(Z);
                        [n,np] = meshgrid(zn,zn);
                        NN = n-np;
                    else
                        Z = zeros(size(en));
                        Z(1) = 1;
                        [n,np] = meshgrid(Z,Z);
                        NN = n-np;
                    end

                    v = squeeze(eigenW(kk,:,:)); % Obtain the corresponding eigen vectors
                    JzhT = JhT.z * ELEf;
                    IzhT = IhT.z * NUCf;
                    tz(:,:,kk,1) = v' * JzhT * v;
                    ttz(:,:,kk,1)  = v' * (JzhT+IzhT) * v;
                    Gc2(:,:,kk,1) = gw0^2 * ttz(:,:,kk,1) .* ttz(:,:,kk,1).' .* NN;
                end
                results.Jmx{iter,iter2,iter3} = tz;
                results.JImx{iter,iter2,iter3} = ttz;
                results.Gc2{iter,iter2,iter3} = Gc2;

                figure;
                hold on
                for ii = 1:N_level
                    plot(fields, sqrt(squeeze(abs(Gc2(ii,ii+1,:)))),'.');
%                     plot(fields, sqrt(abs(squeeze(Gc2(ii+1,ii,:)))),'o'); % symmetrical part of the matrix
                end
                xlabel('Magnetic Field (T)')
                ylabel(ylab)
                title('Nearest level excitation coupling strength')
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
