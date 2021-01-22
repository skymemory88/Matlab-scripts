function [fitresult, gof] = iptopt_0(x, y, field, P0, low, upr)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      Field : current magnetic field
%      P0: Initial guesses for fitting parameters (Ke,w0,G,Br,gamma)
%      Y Output: |S11|
%  Output:
%      fitresult : a fit object representing the S11 parameter.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'S11_inptopt_fit'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( @(kpe, w0, Gc, Br, gma, field, x) mag2db(abs(1+kpe./(1i*(x-w0) - kpe + Gc/(1i*(field-Br)+gma)))),...
    'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'kpe', 'w0', 'Gc', 'Br', 'gma'}, 'problem', {'field'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5)];
opts.Lower = [low(1) low(2) low(3) low(4) low(5)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, 'problem', field, opts);

% Plot fit with data.
figure( 'Name', 'S11_inptopt_fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'S11_inptopt_fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on


