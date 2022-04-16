function [aux_files, location, cmapfig] = fit_list(temp)
if temp == 0.15
    aux_files(1) = "2021_07_0003_0.148K_1.01-4.89T";
%     aux_files(2) = "2021_07_0003_0.142K_7.51-8.49T";
    aux_files(2) = "2021_07_0003_0.148K_6.80-8.79T";
    aux_files(3) = "2021_07_0003_0.142K_8.80-10.19T";
    aux_files(4) = "2021_07_0003_0.148K_10.50-12.19T";
%     aux_files(4) = "2021_07_0003_0.148K_10.31-12.29T";
    aux_files(5) = "2021_07_0003_0.148K_12.20-14.49T";
%     aux_files(5) = "2021_07_0003_0.148K_12.16-13.74T";
    aux_files(6) = "2021_07_0003_0.142K_12.80-17.00T";
%     aux_files(6) = "2021_07_0003_0.142K_13.01-17.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.04\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_150mK_-25dBm_0dB.fig';
elseif temp == 0.18
    aux_files(1) = "2021_07_0004_0.167K_2.00-4.50T";
    aux_files(2) = "2021_07_0004_0.167K_5.51-7.00T";
    aux_files(3) = "2021_07_0004_0.167K_6.95-8.89T";
%     aux_files(3) = "2021_07_0004_0.167K_7.21-9.00T";
    aux_files(4) = "2021_07_0004_0.167K_8.51-10.70T";
    aux_files(5) = "2021_07_0004_0.167K_10.41-12.59T";
    aux_files(6) = "2021_07_0004_0.167K_12.11-13.79T";
%     aux_files(6) = "2021_07_0004_0.167K_12.00-13.79T";
%     aux_files(7) = "2021_07_0004_0.167K_12.50-15.99T";
    aux_files(7) = "2021_07_0004_0.167K_12.71-17.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.05\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_180mK_-25dBm_0dB.fig';
elseif temp == 0.26
    aux_files(1) = "2021_07_0007_0.258K_0.01-4.44T"; % 260 mK
    aux_files(2) = "2021_07_0007_0.258K_5.25-7.13T";
    aux_files(3) = "2021_07_0007_0.258K_7.09-8.94T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.08\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_0-17T_260mK_-25dBm_0dB.fig';
elseif temp == 0.3
    aux_files(1) = "2021_07_0008_0.291K_0.00-4.35T"; % 300 mK
    aux_files(2) = "2021_07_0008_0.291K_5.30-7.33T";
%     aux_files(2) = "2021_07_0008_0.291K_5.51-7.19T";
%     aux_files(2) = "2021_07_0008_0.291K_5.51-7.00T";
    aux_files(3) = "2021_07_0008_0.291K_7.21-9.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.09\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_0-17T_300mK_-25dBm_0dB.fig';
end
end