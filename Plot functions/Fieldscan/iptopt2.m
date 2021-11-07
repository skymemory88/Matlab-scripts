function [fitresult, gof] = iptopt2(x, y, field, P0, low, upr, weight, plt)
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
[xData, yData] = prepareCurveData( x, y );
% Set up fittype and options.
ft = fittype( @(kpe, kpi,  wc, w1, Gc1, w2, Gc2, gma, attn,x) (abs(1+2*kpe./...
    (1i*(x-wc) - (kpe + kpi) + Gc1^2./(1i*(x-w1)-gma) + Gc2^2./(1i*(x-w2)-gma))).*attn),...
    'independent', {'x'}, 'dependent', {'y'},'coefficients', {'kpe', 'kpi', 'wc', 'w1', 'Gc1', 'w2', 'Gc2', 'gma','attn'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)  P0(7)  P0(8) P0(9)];
opts.Lower = [low(1) low(2) low(3) low(4) low(5) low(6) low(7) low(8) low(9)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5) upr(6) upr(7) upr(8) upr(9)];
opts.Weights = weight;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts);

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
