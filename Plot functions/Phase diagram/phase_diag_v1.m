function phase_diag_v1
format long;  
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Phase diagram\functions\');
% addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Phase diagram/functions'));
filepath = 'G:\My Drive\File sharing\PhD program\Research projects\Collaborations\Nb4Rh2C5\Data\2020.09.02';
% filepath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD projects/Collaborations/Daniil Evtushinsky/Transport measurement/2020.03.12';
%The first line is for windows, the second line is for mac OS
%% Choose processing options
opt  = 1;
% file_range = linspace(1,9,9);
file_range = 5;

switch opt
    case 1
        % Simple phase diagram plot from scatter plot
        option1(filepath, file_range)
    case 2
        % Phase diagram plot from interpolation of the raw data + parameter extraction
        option2(filepath, filename)
    case 3
        % Parameter extraction and save files (w/o plots)
        option3(filepath, filename)
end

end

function option1(filepath,file_range)

process = 0;
offset = [0, 0]; % offset to remove instrument artifacts
split = [5.4E-3, 0]; % separate the upper and lower branch

T = cell(size(file_range,2),1); % temperatire
H = cell(size(file_range,2),1); % magnetic field
res = cell(size(file_range,2),1); % raw data
% res_sm = cell(size(file_range,2),1); % filtered data
Tc = cell(size(file_range,2),1); % temperatire
Hc = cell(size(file_range,2),1); % magnetic field
resc = cell(size(file_range,2),1); % raw data
for ii = 1:size(file_range,2)
    filename = sprintf('%04.f',file_range(ii));
%     out = readdata(filepath,filename);
    out = readdataS(filepath,filename);
    
    H{ii} = out.data.Field;
    T{ii} = out.data.Sample_Temperature;
    res{ii} = out.data.Resistance;
%     hall_1(ii,1) = out.data.Lockin1;
%     hall_2(ii,1) = out.data.Lockin2;

    if process == 1
        for jj = 1:length(res{ii})
            if res{ii}(jj) <= 0.0005280
               res{ii}(jj) = res{ii}(jj);
            elseif res{ii}(jj) <= 0.000700
                res{ii}(jj) = res{ii}(jj)-offset(1);
    %             res{ii}(jj) = res{ii}(jj);
            else
                res{ii}(jj) = res{ii}(jj)-offset(2);
            end
        end
    end
%     % Use moving average to smooth the noisy data
%     window = 3;
%     a = 1;
%     b = (1/window)*ones(1,window);
%     res_sm{ii} = filter(b,a,res{ii});
 
    % Define critical field as the mid-point of the transition
    idc = find(round(res{ii},4)==round(mean(res{ii}(res{ii}>=split(1)))/2,4));
    Hc{ii} = mean(H{ii}(idc)); 
    Tc{ii} = mean(T{ii}(idc));
    resc{ii} = mean(res{ii}(idc));
%  
%     figure
%     hold on
%     fig1(ii) = plot(H{ii}(1:5:end),(res{ii}(1:5:end)),'-o'); %Plot field vs. resistance
% %     fig1(ii) = plot(H{ii}(1:300),(res{ii}(1:300)),'-o');
%     plot(Hc{ii},resc{ii},'ro');
%     xlabel ('Nominal magnetic field (T)')
%     ylabel ('Resistance (Ohm)')
%     legend(sprintf('Temp range: %1$.3f - %2$.3f K', min(T{ii}), max(T{ii})));
    
%     figure
%     fig3(ii) = plot(H{ii}(1:10:end),(res_sm{ii}(1:10:end)),'-o');
%     xlabel ('Nominal magnetic field (T)')
%     ylabel ('Smoothed Resistance (Ohm)')
end

% figure
hold on
lgd = strings(1,size(file_range,2));
for ii = 1:size(file_range,2)
    plot(T{ii}(1:2:end),(res{ii}(1:2:end)),'-o'); %Plot temperature vs. resistance
    lgd(ii) = sprintf('Magnetic Field: %.3f T',H{ii}(1));
end
xlabel ('Temperature (K)')
ylabel ('Resistance (Ohm)')
legend(lgd);

temp = cell2mat(T);
field = cell2mat(H);
resistance = cell2mat(res);

Tc = cell2mat(Tc);
Hc = cell2mat(Hc);
resc = cell2mat(resc);

temp_l = min(temp);
temp_h = max(temp);
field_l = min(field);
field_h = max(field);

resistance(isnan(resistance)) = 0; % replace 'NaN' values with 0
    
figure
hold on
scatter3(temp,field,resistance,10,resistance,'filled','MarkerEdgeColor','none');
scatter3(Tc,Hc,resc,10,'r','filled','MarkerEdgeColor','none');
xlabel('Temperature (K)');
ylabel('Magnetic field (T)');
title('2D resistance (raw)');

figure
% plot(Tc,medfilt1(Hc,6),'ro');
plot(Tc,Hc,'ro','MarkerFaceColor','red');
xlabel('Temperature (K)');
ylabel('Magnetic field (T)');
title('Conductance/Resistance phase diagram');

% % Use 2D map to plot interpolated data
% resist = scatteredInterpolant(temp,field,resistance,'nearest','nearest');
% [xq,yq] = meshgrid(linspace(temp_l,temp_h,310),linspace(field_l,field_h,301));
% zq = resist(xq,yq);
% figure
% cmap = pcolor(xq,yq,zq);
% set(cmap, 'edgeColor','none');
% % shading interp;
% colorbar
% xlim([temp_l temp_h])
% ylim([field_l field_h])
% xlabel('Temperature (K)');
% ylabel('Field (T)');

end

function option2(filepath, filename)
end

function option3(filepath, filename)
end