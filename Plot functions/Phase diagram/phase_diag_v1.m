function phase_diag_v1
format long;  
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Phase diagram\functions\');
% addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Phase diagram/functions'));
filepath = 'G:\My Drive\File sharing\PhD program\Research projects\Collaborations\Daniil Evtushinsky\Transport measurement\2020.07.06';
% filepath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD projects/Collaborations/Daniil Evtushinsky/Transport measurement/2020.03.12';
%The first line is for windows, the second line is for mac OS
%% Choose processing options
opt  = 1;
file_range = 10;

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

% T = double.empty(size(file_range,2),321,0); % temperature
% H = double.empty(size(file_range,2),321,0); % field
% hall_1 = double.empty(size(file_range,2),321,0); % Hall measurement 1
% hall_2 = double.empty(size(file_range,2),321,0); % Hall measurement 2
% res = double.empty(size(file_range,2),321,0); % resistrance

T = cell(size(file_range,2),1);
H = cell(size(file_range,2),1);
res = cell(size(file_range,2),1);
for ii = 1:size(file_range,2)
    filename = sprintf('2020_07_%04.f',file_range(ii));
    out = readdata(filepath,filename);

    H{ii} = out.data.DCField1;
    T{ii} = out.data.MC_Temp;
    % hall_1(ii,1) = out.data.Lockin1;
    % hall_2(ii,1) = out.data.Lockin2;
    res{ii} = out.data.Custom;

    figure
    fig(ii) = plot(H{ii}(1:10:end),res{ii}(1:10:end),'-o');
    xlabel ('Nominal magnetic field (T)')
    ylabel ('Hall resistance (Ohm)')
    title('Hall measurement')
end

temp = [T{:}];
field = [H{:}];
resistance = [res{:}];

temp_l = min(temp);
temp_h = max([temp]);
field_l = min([field]);
field_h = max([field]);

resistance(isnan(resistance)) = 0; % replace 'NaN' values with 0
    
figure
hold on
scatter3(temp,field,resistance,10,resistance,'filled','MarkerEdgeColor','none');
xlabel('Temperature (K)');
ylabel('Magnetic field (T)');
title('2D resistivity measurement');

resist = scatteredInterpolant(temp,field,resistance);
[xq,yq] = meshgrid(linspace(temp_l,temp_h,310),linspace(field_l,field_h,301));
zq = resist(xq,yq);

figure
cmap = pcolor(xq,yq,zq);
set(cmap, 'edgeColor','none');
% shading interp;
colorbar
xlim([temp_l temp_h])
ylim([field_l field_h])
xlabel('Temperature (K)');
ylabel('Field (T)');

end

function option2(filepath, filename)
end

function option3(filepath, filename)
end