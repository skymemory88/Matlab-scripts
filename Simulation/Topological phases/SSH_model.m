function [En, rat, N, wvr] = SSH_model(nn, v, plotopt)
% 1D SSH model
N = 2*nn; % force the number of sites to be even
Hr = zeros(N); % real space hamiltonian
sp = length(v)+1; % number of atom species

J = double.empty(0,sp);
J(1,1:end-1) = v; % inter-cell hopping
J(end) = 1.0; % intra-cell hopping

rat = J(1)/J(end); % hopping amplitude ratiod

Options.fluc = 'none'; % hopping amplitude fluctuation
Options.plot = plotopt;
Options.pdc = true;

for ii = 1:sp
    Hr(ii*(N+1) : sp*(N+1) : end) = J(ii);
    Hr(2+(ii-1)*(N+1) : sp*(N+1) : end) = J(ii);
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

if Options.plot == true
    figure;
    hold on
    plot(0:N-1, En,'-s')
    ylabel('Energy')
    yyaxis right
    plot(0:N-1, real(sum(wvr'*wvr)),'-o')
    ylabel('\psi(r)')
    xlabel('Site index')

    figure;
    box on
    hold on
    plot(linspace(0,2*pi,N)-2*pi*floor(linspace(0,2*pi-1e-5,N)./pi), En,'.k');
    plot(-(linspace(0,2*pi,N)-2*pi*floor(linspace(0,2*pi-1e-5,N)./pi)), En,'.k');
    plot([-pi -pi], [-max(En) max(En)],'--r')
    plot([pi pi], [-max(En) max(En)],'--r')
    ylabel('Energy')
    xlabel('Momentum')
end
end