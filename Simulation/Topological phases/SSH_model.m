% SSH model
nn = 200; % half length of the chain
Options.fluc = 'none'; % hopping amplitude fluctuation
Options.plot = true; % Interactive plots
Options.phase = true; % phase diagram
Options.pdc = false; % periodic boundary condition
Options.defc = true; % add a defect to the middle of the chain

% Create a figure object for the interactive plots
fObj1 = figure('WindowState', 'normal', 'Position', [50, 100, 1200, 440]);
ax1 = subplot('Position', [0.1 0.1 0.35 0.75]);
box on
hold(ax1,'on');
[En, N, wvr] = SSHband(nn, 0, Options);
eplot = plot(ax1, 0:N-1, En,'-s');
ylabel('Energy')
yyaxis right
wplot = plot(ax1, 0:N-1, real(sum(conj(wvr)'*wvr)),'--');
ylabel('\psi(r)')
xlabel('Site index')

% Plot the band structure
figure(fObj1);
ax2 = subplot('Position', [0.6 0.1 0.35 0.75]);
box on
hold(ax2, 'on')
plot(ax2, [-pi -pi], [-max(En) max(En)],'--r')
plot(ax2, [pi pi], [-max(En) max(En)],'--r')
rband = plot(ax2, linspace(0,2*pi,N)-2*pi*floor(linspace(0,2*pi-1e-5,N)./pi), En,'.k');
lband = plot(ax2, -(linspace(0,2*pi,N)-2*pi*floor(linspace(0,2*pi-1e-5,N)./pi)), En,'.k');
ylabel('Energy')
xlabel('Momentum')

slider = uicontrol('Parent', fObj1, 'Style', 'slider',...
    'OuterPosition', [100 10 150 20], 'Value', 0, 'min', 0, 'max', 2);
sld.label = uicontrol('Parent',fObj1,'Style','text','Position', [250 10 40 20],...
                'String','J2/J1');

slider.Callback = @(es1,ed) band_update([eplot wplot rband lband], nn, es1.Value, Options);

if Options.phase == true
    phase = figure;
    hold on
    box on
    nn = 12;
    v = 0:0.1:2; % range of inter-cell hoping amplitude
    En = double.empty(2*nn,length(v),0);
    for ii = 1:length(v)
       [En(:,ii,1),~,~] = SSHband(nn, v(ii), Options);
    end
    plot(v,squeeze(En),'-k');
    xlabel('J_2/J_1 ratio')
    ylabel('Energy (arb. unit)')
end

% function to update the interactive plots
function band_update(figs, nn, v, Options)
    [En, ~, wvr] = SSHband(nn, v, Options);
    figs(1).YData = En;
    figs(2).YData = real(sum(conj(wvr)'*wvr));
    figs(3).YData = En;
    figs(4).YData = En;
    drawnow
end

% claculate the bands (nn: half chain length, v: inter-cell bond)
function [En, N, wvr] = SSHband(nn, v, Options)
% 1D SSH model
N = 2*nn; % force the number of sites to be even
Hr = zeros(N); % real space hamiltonian
sp = length(v)+1; % number of atom species

J = double.empty(0,sp);
J(1,1:end-1) = v; % inter-cell hopping
J(end) = 1.0; % intra-cell hopping

for ii = 1:sp
    Hr(ii*(N+1) : sp*(N+1) : end) = J(ii);
    Hr(2+(ii-1)*(N+1) : sp*(N+1) : end) = J(ii);
end

if Options.defc == true
    mid = round(length(Hr)/2);
    utemp = Hr(mid,mid+1);
    ltemp = Hr(mid+1,mid);
    Hr(mid, mid+1) = Hr(mid+1,mid+2);
    Hr(mid+1, mid) = Hr(mid+2,mid+1);
    Hr(mid+1,mid+2) = utemp;
    Hr(mid+2,mid+1) = ltemp;
end

switch Options.fluc
    case 'random'
        rng(0,'twister'); % seed the random number generator
        bumps = J2 + rand(1,nn-1);
        Hr(2*(N+1) : sp*(N+1) : end) = J(2) + bumps;
        Hr(2+(ii-1)*(N+1) : sp*(N+1) : end) = J(2) + bumps;
    case 'uniform'
        rng(0,'twister'); % seed the random number generator
        bumps = rand(1)/10;
        Hr(2*(N+1) : sp*(N+1) : end) = J(2) + bumps;
        Hr(2+(ii-1)*(N+1) : sp*(N+1) : end) = J2 + bumps;
end
% periodic boundary conditiond
if Options.pdc == true
    Hr(1,N) = J(end);
    Hr(N,1) = J(end);
end
[wvr, Er] = eig(Hr);
En = diag(Er);
end