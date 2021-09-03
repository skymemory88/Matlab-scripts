function Perturb
% close all
hold on
clearvars
temp = 0.130;
kB = 1/11.6; %Boltzmann constant [meV/K]
beta = 1/(temp*kB);

Options.Elevel = false;
Options.Ediff = true;
Options.savedata = false;

field_l = 0;
field_h = 9;
H0 = linspace(field_l,field_h,300); % sampling points along field axis
w0 = 3.675; % Cavity resonant freqeuency
g = 0.001; % Coupling strength measured against the ground state energy
% g = 0.01*w0;

% color = ["black","red","blue","magenta","green","yellow","cyan"];
theta = 0; % Angle (in degrees) in the transverse field plane
phi = 0;
% marker = [":","-.","--","-"];
% lg = strings(1,2); % Create an empty array for legends
% figs = gobjects(7,numel(temp));

for iter = 1:numel(temp)
    for iter2 = 1:numel(theta)
%         cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Mean Field --Peter & Ivan\output\A=1.0_angles_Yikai');
        addpath(genpath('G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations results\Matlab\Susceptibilities\without Hz_I'));
%         cd('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Simulation/External source/Mean Field --Peter & Ivan/output/A=1.0_angles_Peter')
        lname=['Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$.2fDeg_%3$.1fDeg',temp(iter),theta,phi),'.mat'];
        
        load(lname,'-mat','eee','fff');
        fields = vecnorm(fff);
        E2f = 241.8; % Convert Energy to frequency
        E(:,:)=squeeze(eee)*E2f;
        %E(:,:) = eee(:,1,:)*E2F;
        Ediff = double.empty(7,size(E,1),0);       
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

            bzF = double.empty(size(Ediff,1),size(fields,2),0);
            for ii=1:7 % To use all 7 transition lines
                Ediff(ii,:,1)=E(:,ii+1)-E(:,ii);
                bzF(ii,:,1) = exp(-Ediff(ii,:)*beta); % Boltzmann factor for each transition line
            end
%             figs(:, iter2) = plot(fields, Ediff, 'Marker', 'none', 'LineStyle',marker(1), 'Color',color(iter2),'LineWidth',2);
%             hold on
%             lg(iter, iter2) = [num2str(temp(iter)*1000,"T = %u mK, A = A_{th}"),num2str(theta(iter2),", theta = %u")];
%             set(gca,'fontsize',15,'Xtick',0:1:max(fields));
%             xlabel('Field (T)','FontSize',15);
%             ylabel('Energy difference (GHz)','FontSize',15);
%             grid off
%             box on;
%             xlim([min(fields) max(fields)]);
%             ylim([0.9*min(Ediff,[],'All') 1.1*max(Ediff,[],'All')]);
%             title('Difference between energy levels','FontSize',15);
        end
        %% Interpolate the dispersion relations
%         % Test dispersion relations (7 straight lines)
%         a = -0.2;
%         b = linspace(4.2,4.9,8)';
%         w1 = repmat(a*H0,8,1);
%         w = w1 + b.*ones(8,size(H0,2));
%         w(1,:) = w0;
        
        % Calcualted dispersion relations
        w = double.empty(size(Ediff,1)+1,size(H0,2),0);
        boltFac = double.empty(size(Ediff,1),size(H0,2),0);
        w(1,:,1) = w0;
        for ii = 2:size(Ediff,1)+1
            w(ii,:) = interp1(fields,Ediff(ii-1,:),H0); 
            boltFac(ii-1,:,1) = interp1(fields,bzF(ii-1,:),H0);
        end
        ww = double.empty(0,size(Ediff,1)+1,size(H0,2));
%         for jj = 1:size(H0,2) % trial Hamiltonian 1
%             gs = g*abs(w0-w(2:end,jj));
% %             gs = g*abs(w0-w(2,jj))*ones(1,size(Ediff,1));
%             gs(1:2) = 0;
%             gs(4:end) = 0;
%             Ham = diag(gs,1);
%             Ham = Ham + Ham';
%             Ham = Ham + w(:,jj).*eye(size(w,1));
%             ww(1,:,jj) = eig(Ham);
%         end

        for jj = 1:size(H0,2) % Trial Hamiltonian 2 (Correct one!)
            gs = g./abs(w0-w(2:end,jj));
%             gs = g/abs(w0-w(2,jj))*ones(1,size(Ediff,1))*(boltFac(:,jj)/max(boltFac(:,jj)));
%             gs = g.*(boltFac(:,jj)/max(boltFac(:,jj)));
%             gs(1:2) = 0;
%             gs(1:end-1) = 0;
            Ham = diag(w(:,jj));
            Ham(1,2:end) = gs;
            Ham(2:end,1) = gs';
            ww(1,:,jj) = eig(Ham);
        end
        
        ww = squeeze(ww);
%         figure
        pert = plot(H0,ww,'--k','linewidth',1);
        xlim([field_l field_h])
        ylim([3.65 3.70])
        delete(pert)
    end
end