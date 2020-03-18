clear
data = importdata('table_0K_0T.txt','\t',1);
% Data import, column index should be adjusted base on the raw data
num = data.data(:,5);
E_dens = data.data(:,7);
temp = data.data(:,12);
Bz = data.data(:,10);
Q = data.data(:,11);

%% Option 1: Plot the energy density versus Nx from a single scan with the minimum
%marked
[E_min, loc] = findpeaks(-E_dens);
plot(num/1e9,E_dens,'o-')
hold on
plot(num(loc)/1e9,E_dens(loc),'ko','MarkerFaceColor','r')
xlabel('Skyrmion lattice constant (nm)');
ylabel('Energy density (J/m^3)');

%% Option 2: Plot the phase diagram from interpolations from raw data
% %make a grid for contour plot (5 times of the density of the raw data)
% xl = linspace(min(temp),max(temp),5*((max(temp)-min(temp))/TStep)); 
% yl = linspace(min(Bz),max(Bz),5*((max(Bz)-min(Bz))/BStep)) ; 
% [X,Y] = meshgrid(xl,yl);
% 
% %Interpolate data on the mesh grid
% coords = [temp,Bz];
% intrp = scatteredInterpolant(coords,Q);
% intrp.Method = 'nearest';
% Q_intrp = intrp([X(:),Y(:)]);
% Q_intrp = reshape(Q_intrp,size(X));
% 
% %Plot the contour map
% contourf(X,Y,Q_intrp);
% xlabel('Temperature (K)');
% ylabel('Magnetic field (T)');
% colorbar('eastoutside');
% 
% %Plot the raw data for comparision
% figure
% scatter3(temp,Bz,Q);
% xlabel('Temperature (K)');
% ylabel('Magnetic field (T)');

% % Find the local minimum in the energy density curve as a function of
% % lattice constant with manual set Nx ranges, used for multi-scan data
% % WARNING: may generate a large number of figures !!
% NxMin = 70;
% NxStep = 2;
% NxMax = 188;
% NxRange = (NxMax-NxMin)/NxStep;
%
% for lineNum = 1:NxRange+1:length(num)-NxRange
%     [E_min, loc] = findpeaks(-E_dens(lineNum:lineNum+NxRange));
%     figure
%     plot(num(lineNum:lineNum+NxRange)/1e-9,E_dens((lineNum:lineNum+NxRange)),'o-')
%     hold on
%     plot(num(loc)/1e-9,-E_min,'ko','MarkerFaceColor','r')
%     hold off
%     xlabel('Skyrmion lattice constant (nm)');
%     ylabel('Energy density (J/m^3)');
% end

% % Exract phase from minimum energy calculation and mark down the phase as
% "+1"
% for lineNum = 1:NxRange+1:length(num)-NxRange
%     [E_min, loc] = findpeaks(-E_dens(lineNum:lineNum+NxRange));
%     density(lineNum + loc) = density(lineNum + loc) + 1;
%     if loc ~= 0
%         density(lineNum) = density(lineNum) + 1;
%     end
% end
% 
% density = zeros(mod(length(num),NxRange));
% phase = density(1:NxRange+1:length(num)-NxRange);
% Bz2 = Bz(1:NxRange+1:length(num)-NxRange);
% plot(Bz2,phase,'o-');
% xlabel('Mangetic Field (T)')
% ylabel('Skyrminion phase')
% legend('Skyrmion lattice at 0 K')