% susceptibility plots
clearvars freq_total fields chir chii chi idx

temp = 0.150; % temperature
theta = 0.0; % field angle in ac-plane
phi = 2.3; % field angle in ab-plane
gama = 9e-5; % spin line width

Option.hyp = true; % hyperfine interaction option
Option.RPA = false; % RPA option
Option.nZee = true; % nuclear Zeeman interaction option
Option.sPlot = true; % isolated plots of freqeuency slice
    part = 'both'; % Plot the 'real', 'imaginary', or 'both' parts of susceptibility
    pFormat = 'condensed'; % 1.'condensed': on one plot; 2. 'separate': individual plots; 3. 'none'
Option.peaks = true; % extract the peak positions of Im[X]

B0 = [5.015]; % field location of frequency cuts to plot/fit

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

file_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\',...
    'Research projects\LiHoF4 project\Data\Simulations\Matlab\Susceptibilities'];
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
ReSlice = gobjects(length(nonzeros(chip)), 1);
ImSlice = gobjects(length(nonzeros(chip)), 1);

idx = 1;
for ii = 1:length(chip(:))
    if chip(ii) ~= 0
        clearvars sz fsz lsz cb
        [row, col] = ind2sub(size(chip),ii);
        switch part
            case 'real'
                chir = real(squeeze(chi(row,col,:,:)));
            case 'imaginary'
                chii = imag(squeeze(chi(row,col,:,:)));
            case 'both'
                chir = real(squeeze(chi(row,col,:,:)));
                chii = imag(squeeze(chi(row,col,:,:)));
        end
        
        bidx = uint32.empty(length(B0),0); % container for field indices
        for jj = 1:length(B0)
            [~, bidx(jj)] = min(abs(B0(jj)-fields(:))); % locate the data slices
        end
        
        widx = uint32.empty(length(fields),7,0); % container for frequency indices
        w0 = double.empty(length(fields),7,0); % container for frequency indices
        for kk = 1:length(fields)
           [~, fidx] = findpeaks(abs(chii(:,kk)), 'NPeaks', 7);
           if isempty(fidx)
               w0(kk,:,1) = 0;
               widx(kk,:,1) = NaN;
           else
               f0 = freq_total(fidx);
               f0(end+1:7) = 0;
               widx = fidx;
               w0(kk,:,1) = f0;
           end
        end
        
        if exist('chir', 'var') && Option.sPlot == true
            figure
            cmp = pcolor(fields, freq_total, chir);
            cmp.EdgeColor = 'none';
            fsz = get(gca, 'Position');
            colormap(cmap)
%             caxis([-40 120])
            xlabel('Magnetic Field (T)');
            ylabel('Frequency (GHz)')
            set(gca,'FontSize',14)
            xticks([0:2:17])
            
%             colorbar off
%             title ''
            title(strcat("Re[\chi_", tsr(row,col),']'))
            cbr = colorbar('Location','Northoutside');
            sz = get(cbr, 'Position');
            set(cbr,'Position',[sz(1) sz(2) sz(3)-0.15 sz(4)]);
            set(gca, 'Position', [fsz(1) fsz(2) fsz(3) fsz(4)-0.15]);
            
            cbr.Label.String = 'Amp.';
            lsz = cbr.Label.Position;
%             cb.Label.Position = [lsz(1)*2.2 lsz(2)-2.25 0];
            
            for jj = 1:length(B0)
                switch pFormat
                    case 'condensed'
                        if idx == 1
                            ReSlice(idx) = figure;
                        else
                            figure(ReSlice(1));
                        end
                    case 'separate'
                        ReSlice(idx) = figure;
                    case 'none'
                        break
                end
                box on
                hold on
                plot(freq_total, chir(:,bidx(jj)), 'LineStyle', '-', 'LineWidth', 2);
                xlabel('Frequency (GHz)')
                ylabel(strcat("Re[\chi_", tsr(row,col),"]"))
                set(gca, 'FontSize', 14)
            end
        end
        
        if exist('chii', 'var') && Option.sPlot == true
            figure
            cmp = pcolor(fields, freq_total, chii);
            cmp.EdgeColor = 'none';
            fsz = get(gca, 'Position');
            colormap(cmap)
%             caxis([-40 120])
            xlabel('Magnetic Field (T)');
            ylabel('Frequency (GHz)')
            set(gca,'FontSize',14)
            xticks([0:2:17])
            
%             colorbar off
%             title ''
            title(strcat("Im[\chi_", tsr(row,col),']'))
            cbi = colorbar('Location','Northoutside');
            sz = get(cbi, 'Position');
            set(cbi,'Position',[sz(1) sz(2) sz(3)-0.15 sz(4)]);
            set(gca, 'Position', [fsz(1) fsz(2) fsz(3) fsz(4)-0.15]);
            
            cbi.Label.String = 'Amp.';
            lsz = cbi.Label.Position;
%             cb.Label.Position = [lsz(1)*2.2 lsz(2)-2.25 0];
            
            for jj = 1:length(B0)
                switch pFormat
                    case 'condensed'
                        if idx == 1
                            ImSlice(idx) = figure;
                            box on
                            hold on
                        else
                            figure(ImSlice(1));
                        end
                    case 'separate'
                        ImSlice(idx) = figure;
                    case 'none'
                        break
                end
                plot(freq_total, chii(:,bidx(jj)), 'LineStyle', '-', 'LineWidth', 2);
                xlabel('Frequency (GHz)')
                ylabel(strcat("Im[\chi_", tsr(row,col),"]"))
                set(gca, 'FontSize', 14)
            end
        end
        idx = idx + 1;
    end
end