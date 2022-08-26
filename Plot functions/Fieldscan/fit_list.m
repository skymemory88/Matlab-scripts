function [aux_files, location, cmapfig] = fit_list(temp)
if temp == 0.10
    aux_files(1) = "2021_07_0003_0.148K_1.01-4.89T";
    aux_files(2) = "";
    aux_files(3) = "2021_07_0003_0.148K_6.80-8.79T";
%     aux_files(3) = "2021_07_0003_0.142K_7.51-8.49T"; % alternative/additional fits
    aux_files(4) = "2021_07_0003_0.142K_8.80-10.19T";
    aux_files(5) = "2021_07_0003_0.148K_10.50-12.19T";
%     aux_files(5) = "2021_07_0003_0.148K_10.31-12.29T"; % alternative/additional fits
    aux_files(6) = "2021_07_0003_0.148K_12.20-13.69T";
%     aux_files(6) = "2021_07_0003_0.148K_12.20-15.00T"; % alternative/additional fits
%     aux_files(6) = "2021_07_0003_0.148K_12.20-14.49T"; % alternative/additional fits
%     aux_files(7) = "2021_07_0003_0.142K_12.80-17.00T";
    aux_files(7) = "2021_07_0003_0.142K_13.01-17.00T"; % alternative/additional fits
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.04\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_150mK_-25dBm_0dB.fig';
elseif temp == 0.135
    aux_files(1) = "2021_07_0004_0.167K_2.00-4.50T";
    aux_files(2) = "2021_07_0004_0.167K_5.40-7.15T";
%     aux_files(2) = "2021_07_0004_0.167K_5.31-7.09T"; % alternative/additional fits
%     aux_files(2) = "2021_07_0004_0.167K_5.51-7.00T"; % alternative/additional fits
    aux_files(3) = "2021_07_0004_0.167K_6.95-8.89T";
%     aux_files(3) = "2021_07_0004_0.167K_7.10-9.00T"; % alternative/additional fits
    aux_files(4) = "2021_07_0004_0.167K_8.51-10.70T";
    aux_files(5) = "2021_07_0004_0.167K_10.21-12.09T";
%     aux_files(5) = "2021_07_0004_0.167K_10.41-12.59T"; % alternative/additional fits
    aux_files(6) = "2021_07_0004_0.167K_12.11-13.79T";
%     aux_files(6) = "2021_07_0004_0.167K_12.41-14.49T"; % alternative/additional fits
%     aux_files(6) = "2021_07_0004_0.167K_12.00-13.79T"; % alternative/additional fits
    aux_files(7) = "2021_07_0004_0.167K_12.71-17.00T";
%     aux_files(7) = "2021_07_0004_0.167K_12.50-15.99T"; % alternative/additional fits
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.05\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_17-0T_180mK_-25dBm_0dB.fig';
elseif temp == 0.18
    aux_files(1) = "2021_07_0006_0.201K_0.00-4.71T";
    aux_files(2) = "2021_07_0006_0.201K_5.29-7.16T";
%     aux_files(2) = "2021_07_0006_0.201K_5.25-7.19T"; % alternative/additional fits
    aux_files(3) = "2021_07_0006_0.201K_7.21-8.99T";
%     aux_files(3) = "2021_07_0006_0.201K_7.08-8.89T"; % alternative/additional fits
    aux_files(4) = "2021_07_0006_0.201K_8.60-10.80T";
%     aux_files(4) = "2021_07_0006_0.201K_8.81-10.80T"; % alternative/additional fits
    aux_files(5) = ""; % space filler
    aux_files(6) = ""; % space filler
    aux_files(7) = "2021_07_0006_0.201K_12.71-15.99T";
%     aux_files(7) = "2021_07_0006_0.201K_13.01-15.99T"; % alternative/additional fits
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.07\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_0-17T_220mK_-25dBm_0dB.fig';
elseif temp == 0.24
    aux_files(1) = "2021_07_0007_0.258K_0.01-4.44T";
    aux_files(2) = "2021_07_0007_0.258K_5.51-7.19T";
%     aux_files(2) = "2021_07_0007_0.258K_5.25-7.13T"; % alternative/additional fits
    aux_files(3) = "2021_07_0007_0.258K_7.25-8.91T";
%     aux_files(3) = "2021_07_0007_0.258K_7.09-8.94T"; % alternative/additional fits
    aux_files(4) = ""; % space filler
    aux_files(5) = ""; % space filler
    aux_files(6) = ""; % space filler
    aux_files(7) = "2021_07_0007_0.258K_12.61-15.99T";
%     aux_files(7) = "2021_07_0007_0.258K_12.17-14.70T"; % alternative/additional fits
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.08\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_0-17T_260mK_-25dBm_0dB.fig';
elseif temp == 0.25
    aux_files(1) = "2021_07_0008_0.291K_0.00-4.35T"; % 290 mK
    aux_files(2) = "2021_07_0008_0.291K_5.31-7.09T";
%     aux_files(2) = "2021_07_0008_0.291K_5.31-7.39T"; % alternative/additional fits
%     aux_files(2) = "2021_07_0008_0.291K_5.30-7.33T"; % alternative/additional fits
%     aux_files(2) = "2021_07_0008_0.291K_5.51-7.19T"; % alternative/additional fits
%     aux_files(2) = "2021_07_0008_0.291K_5.51-7.00T"; % alternative/additional fits
    aux_files(3) = "2021_07_0008_0.291K_7.21-9.00T";
    aux_files(4) = ""; % space filler
    aux_files(5) = ""; % space filler
    aux_files(6) = ""; % space filler
    aux_files(7) = "2021_07_0008_0.291K_13.01-15.00T";
    location = "G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC239\2021.07.09\";
    cmapfig = 'LiHoF4_3.52-3.76GHz_0-17T_300mK_-25dBm_0dB.fig';
end
end