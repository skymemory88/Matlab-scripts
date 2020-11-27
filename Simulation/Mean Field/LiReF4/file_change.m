function file_change()
Temperatures = [0.1 0.12 0.15 0.2 0.24 0.3 0.35 0.4 0.45 0.5 0.8 1.2 1.6 1.8 2];
% Temperatures = 0.08;
phi = 0; % ab-plane rotation, phi=0 means H along x (in radian)
cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output\without Hz_I');
    for ii = 1:length(Temperatures)
        filenames = strcat('Hscan_LiHoF4_', sprintf('%1$3.3fK_%2$uDeg', Temperatures(ii), phi),'.mat');
        load(filenames,'-mat','field');
        fff = field;
        save(filenames,'-mat','fff','-append');
        disp([filenames, ' changed']);
    end
end