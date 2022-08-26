% susceptibility plots
clearvars chir chii chi

temp = 0.3; % temperature
theta = 0.0; % field angle in ac-plane
phi = 0.0; % field angle in ab-plane
gama = 9e-5; % spin line width

Option.hyp = true; % hyperfine interaction option
Option.RPA = true; % RPA option
Option.nZee = false; % nuclear Zeeman interaction option
Option.sPlot = true; % isolated plots of freqeuency slice
Option.part = 'both'; % choose the 'real', 'imaginary', or 'both' part to plot

B0 = [3.825]; % field location of frequency cuts to plot/fit
bidx = uint32.empty(length(B0),0);

cmap = unique([[0 0 0];[zeros(50,1),linspace(0,0.5,50)',linspace(0.4,1,50)'];...
               [linspace(0,1,70)',0.5*ones(70,1),linspace(1,0,70)'];...
               [ones(60,1),linspace(0.5,1,60)',linspace(0,1,60)']],'rows');

mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250];...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool

chip = [0 0 0
        0 0 0
        0 0 1]; % choose the susceptibility tensor element to plot

tsr = ["{xx}"  "{xy}"  "{xz}"
       "{yx}"  "{yy}"  "{yz}"
       "{zx}"  "{zy}"  "{zz}"];

% file_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\',...
%     'Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities'];
file_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Simulations\',...
    'Matlab\Susceptibilities\'];
file_part = sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_%4$.2e_*GHz', temp, theta, phi, gama);

if Option.nZee == true
    nZeePath = '\with Hz_I\';
else
    nZeePath = '\without Hz_I\test\';
end

if Option.RPA ==true
    filename = strcat('RPA_LiHoF4_', file_part, ".mat");
else
    filename = strcat('chi0_LiHoF4_', file_part, ".mat");
end

file_loc = [file_loc, nZeePath];
chi_file = dir( fullfile(file_loc, filename) );
load([file_loc chi_file.name],'-mat','fields','freq_total','chi','unit'); % load calculated suscpetibilities
slice = gobjects(length(nonzeros(chip)), 1);

idx = 1;
for ii = 1:length(chip(:))
    if chip(ii) ~= 0
        clearvars sz fsz lsz cb

        [row, col] = ind2sub(size(chip),ii);
        switch Option.part
            case 'real'
                chir = real(squeeze(chi(row,col,:,:)));
            case 'imaginary'
                chii = imag(squeeze(chi(row,col,:,:)));
            case 'both'
                chir = real(squeeze(chi(row,col,:,:)));
                chii = imag(squeeze(chi(row,col,:,:)));
        end

        if exist('chir', 'var')
            figure
            cmp = pcolor(fields, freq_total, chir);
            cmp.EdgeColor = 'none';
            fsz = get(gca, 'Position');
            colormap(cmap)
    %         caxis([-40 120])
            xlabel('Magnetic Field (T)');
            ylabel('Frequency (GHz)')
            set(gca,'FontSize',14)
            xticks([0:2:17])
    
    %         colorbar off
    %         title ''
            title(strcat("\chi_", tsr(row,col)))
            cb = colorbar('Location','Northoutside');
            sz = get(cb, 'Position');
            set(cb,'Position',[sz(1) sz(2) sz(3)-0.15 sz(4)]);
            set(gca, 'Position', [fsz(1) fsz(2) fsz(3) fsz(4)-0.15]);
    
            cb.Label.String = 'Amp.';
    %         lsz = cb.Label.Position;
    %         cb.Label.Position = [lsz(1)*2.7 lsz(2)-2.25 0];
            if Option.sPlot == true
                for jj = 1:length(B0)
                    [~, bidx(jj)] = min(abs(B0(jj)-fields(:)));
                    slice(idx) = figure;
                    box on
                    hold on
                    plot(freq_total, chir(:,bidx(jj)), 'LineStyle', '-');
                end
                xlabel('Frequency (GHz)')
                ylabel(strcat("Re[\chi_", tsr(row,col),"]"))
            end
        end
        
        if exist('chii', 'var')
            figure
            cmp = pcolor(fields, freq_total, chii);
            cmp.EdgeColor = 'none';
            fsz = get(gca, 'Position');
            colormap(cmap)
    %         caxis([-40 120])
            xlabel('Magnetic Field (T)');
            ylabel('Frequency (GHz)')
            set(gca,'FontSize',14)
            xticks([0:2:17])
    
    %         colorbar off
    %         title ''
            title(strcat("\chi_", tsr(row,col)))
            cb = colorbar('Location','Northoutside');
            sz = get(cb, 'Position');
            set(cb,'Position',[sz(1) sz(2) sz(3)-0.15 sz(4)]);
            set(gca, 'Position', [fsz(1) fsz(2) fsz(3) fsz(4)-0.15]);
    
            cb.Label.String = 'Amp.';
    %         lsz = cb.Label.Position;
    %         cb.Label.Position = [lsz(1)*2.7 lsz(2)-2.25 0];

            if Option.sPlot == true
                for jj = 1:length(B0)
                    [~, bidx(jj)] = min(abs(B0(jj)-fields(:)));
                    slice(idx) = figure;
                    box on
                    hold on
                    plot(freq_total, chii(:,bidx(jj)), 'LineStyle', '-');
                end
                xlabel('Frequency (GHz)')
                ylabel(strcat("Im[\chi_", tsr(row,col),"]"))
            end
        end

        idx = idx + 1;
    end
end