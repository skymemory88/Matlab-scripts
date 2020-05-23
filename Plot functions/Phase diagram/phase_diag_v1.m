function phase_diag_v1
format long;  
addpath('G:\My Drive\File sharing\Programming scripts\Matlab\Plot functions\Phase diagram\functions\');
% addpath(genpath('/Volumes/GoogleDrive/My Drive/File sharing/Programming scripts/Matlab/Plot functions/Phase diagram/functions'));
filepath = 'G:\My Drive\File sharing\PhD projects\Collaborations\Daniil Evtushinsky\Transport measurement\2020.03.12';
% filepath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD projects/Collaborations/Daniil Evtushinsky/Transport measurement/2020.03.12';
%The first line is for windows, the second line is for mac OS
filename = '2020_03_0004';

%% Choose processing options
opt  = 1;

switch opt
    case 1
        % Simple phase diagram plot from scatter plot
        option1(filepath, filename)
    case 2
        % Phase diagram plot from interpolation of the raw data + parameter extraction
        option2(filepath, filename)
    case 3
        % Parameter extraction and save files (w/o plots)
        option3(filepath, filename)
end

end

function option1(filepath,filename)
out = readdata(filepath,filename);

H = out.data.DCField1;
T = out.data.Sample_Temp;
lck = out.data.Lockin1;
res = out.data.Resistivity;

h_low = find(H==min(H));
h_high = find(H==max(H));

figure
hold on
plot(H(h_low(1):h_high(1)),lck(h_low(1):h_high(1)),'o');
xlabel ('Nominal magnetic field (T)')
ylabel ('Hall resistance (Ohm)')
title('Hall measurement')

res(isnan(res)) = 0; % replace 'NaN' values with 0

figure
hold on
scatter3(T,H,res,10,res,'filled','MarkerEdgeColor','none');
xlabel('Temperature (K)');
ylabel('Magnetic field (T)');
title('2D resistivity measurement');

end

function option2(filepath, filename)
end

function option3(filepath, filename)
end