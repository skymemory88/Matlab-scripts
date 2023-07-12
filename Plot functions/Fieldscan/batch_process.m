% Operation options
% 1. Simple color plot of minimally processed raw data
% 2. Color plot + Data fitting
% 3. Data fitting and file saving (w/o plots)

opt.type = 1;
opt.nZVL = 1;
opt.resnt = false;

% Figure plot options:
plotopt.cplfit = false;
plotopt.lnwd = 2;
plotopt.ftsz = 12;
plotopt.mksz = 2;

% File locations
filepath = 'G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC200\2020.08.25';
%     filepath = '/Volumes/GoogleDrive/My Drive/File sharing/PhD program/Research projects/LiHoF4 project/Data/Experiment/LiHoF4/SC127/SC127_2 (2.5 x 1 x 0.5 mm, triangle)/25.07.2020';
%The first line is for windows, the second line is for mac OS

file_range = 44;
for ii = 1:length(file_range)
    filename = sprintf('2020_08_00%u.dat',file_range(ii));
    fileobj = fullfile(filepath,filename);
    Field_scan_v5a(fileobj,opt,plotopt);
end

clearvars

