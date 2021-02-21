function [fitresult, gof] = iptopt_0(x, y, field, B0, spin, P0, low, upr, weight, plt)
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
%% Fit: 'S11_inptopt_fit'.
[xData, yData] = prepareCurveData( x, y );
% Set up fittype and options.
ft = fittype( @(kpe, kpi, w0, Gc, gma, shift, field, B0, spin, x) (shift + abs(1+kpe./...
    (1i*(x-w0) - (kpe + kpi) + Gc^2./(1i*(x-(spin*(field-B0)+w0) )-gma)))),...
    'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'kpe', 'kpi', 'w0', 'Gc', 'gma', 'shift'}, 'problem', {'field','B0', 'spin'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)];
opts.Lower = [low(1) low(2) low(3) low(4) low(5) low(6)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5) upr(6)];
opts.Weights = weight;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'problem', {field, B0, spin}, opts);

if plt
    % Plot fit with data.
    figure( 'Name', sprintf('S11_inptopt_fit at B = %.2f T',field) );
    h = plot( fitresult, xData, yData );
    legend( h, 'Data', 'S11_fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'Frequency (GHz)', 'Interpreter', 'none' );
    ylabel( '|S11|', 'Interpreter', 'none' );
    grid on
end
