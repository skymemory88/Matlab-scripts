function [fitresult, gof] = iptopt(x, y, coef, P0, low, upr, weight, plt)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      alpha : coefficient of linearized w0(B)
%      P0: Initial guesses for fitting parameters (Ke,Ki,w0,Gc,gamma)
%      Y Output: |S11|
%      Plt (true of false): Plot the fit along with data
%  Output:
%      fitresult : a fit object representing the S11 parameter.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
% Fit: 'S11_iptopt_fit'.
mu0 = 4*pi*10^-7; % Vacuum permeability ([H/m])
hbar = 1.055E-34; % Reduced Planck constant [J.s]
meV2J = 1.602217e-22; % [J/meV]
rho = 4e30/(5.17*5.17*10.75); % Holmium (magnetic moment) number density [m^-3]
f2E = hbar*2*pi*10^9/meV2J; % [meV/GHz]

kpe = coef(1); % external dissipation rate
kpi = coef(2); % internal dissipation rate
wc = coef(3); % cavity resonance frequency
field = coef(4); % current magnetic field
fill = coef(5); % Filling factor (Sample specific)

[xData, yData] = prepareCurveData( x, y );
g = sqrt(mu0*wc*10^9*2*pi*rho*fill/hbar/2); % susceptibility prefactor [T.(J.s)^-1]
g = g * meV2J * f2E * 10^-9;
% Set up fittype and options.
ft = fittype( @(xr_off, xi_off, omega, Gc, gma, attn, kpe, kpi, wc, x) (abs(1 + 2*kpe./...
    (1i*(x - wc) - (kpi + kpe) + 1i*g^2*(xr_off + 1i*xi_off + Gc./(1i*(x-omega )-gma)))).*attn),...
    'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'xr_off', 'xi_off', 'omega', 'Gc', 'gma','attn'},...
    'problem',{'kpe','kpi','wc'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)];
opts.Lower = [low(1) low(2) low(3) low(4) low(5) low(6)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5) upr(6)];
opts.Weights = weight;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc});

if plt
    % Plot fit with data.
    figure( 'Name', sprintf('S11_inptopt_fit at B = %.3f T',field) );
    h = plot( fitresult, xData, yData );
    legend( h, 'Data', 'S11_fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'Frequency (GHz)', 'Interpreter', 'none' );
    ylabel( '|S11|', 'Interpreter', 'none' );
    grid on
end
