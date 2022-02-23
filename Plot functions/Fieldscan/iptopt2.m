function [fitresult, gof] = iptopt2(x, y, coef, param, aux_param, weight, plt)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      P0: Initial guesses for fitting parameters (Ke,Ki,wc,Gc,gamma)
%      Y Output: |S11|
%      Plt (true of false): Plot the fit along with data
%  Output:
%      fitresult : a fit object representing the S11 parameter.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'S11_inptopt_fit'.
kpe = coef(1); % external dissipation rate
kpi = coef(2); % internal dissipation rate
wc = coef(3); % cavity resonance frequency
field = coef(4); % current magnetic field

w2 = aux_param.f0(2); % auxilliary spin resonance frequency
% gc0 = aux_param.gc(2); % auxilliary spin resonance coupling strength

P0 = param(1,:); % starting values of the fitting parameters
low = param(2,:); % lower limit of the fitting parameters
upr = param(3,:); % upper limit of the fitting parameters

[xData, yData] = prepareCurveData( x, y );
% Set up fittype and options.
ft = fittype( @(omega, Gc, gma, attn, xr, xi, Gc2, kpe, kpi, wc, w2, x) (abs(1+2*kpe./...
    (1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma) + Gc2^2./(1i*(x-w2)-gma))).*attn),...
    'independent', {'x'}, 'dependent', {'y'},'coefficients', {'omega', 'Gc',  'gma','attn', 'xr', 'xi', 'Gc2'},...
    'problem', {'kpe','kpi','wc','w2'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1)   P0(2)   P0(3)   P0(4)   0    0   P0(2)]; % {'omega', 'Gc', 'gma', 'attn', 'xr', 'xi', 'Gc2'}
opts.Lower = [low(1)   low(2)  low(3)  low(4)    0    0  low(2)];
opts.Upper = [upr(1)   upr(2)  upr(3)  upr(4)    0    0  upr(2)];
opts.Weights = weight;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {kpe, kpi, wc, w2});

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
