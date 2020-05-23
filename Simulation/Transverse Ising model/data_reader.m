mpwd = 13.9e-3;

% molar mass
molarmass = 245.43;

% number of magnetic ions/formula
Nions = 1;

% number of mols of magnetic ions
Nmols = Nions*(mpwd/molarmass);

% measuring field
Hmsr = 100;    % in Oe

LiDyF4 = importdata('G:\My Drive\Programming scripts\Matlab/0066.dc.dat',',',31);

xpwd = LiDyF4.data(:,4);
ypwd = LiDyF4.data(:,5)/Hmsr/Nmols; % in emu/mol Oe

plot(xpwd,ypwd)
%plot(xpwd,1./ypwd)