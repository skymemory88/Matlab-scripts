function [fitresult, gof] = wk_cpl_fit(H0, f0, spin, hPara, fPara)
%CREATEFIT(H0,F0)
%  Create a fit.
%
%  Data for 'wk_cpl_fit' fit:
%      X Input : H0
%      Y Output: f0
%      hPara: [init_h, hBound_l, hBound_h]
%      fPara: [init_f, fBound_l, fBound_h]
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'wk_cpl_fit'.
[xData, yData] = prepareCurveData( H0, f0 );

% Set up fittype and options.
ft = fittype( 'wc - g^2*spin*(x-x0)/((spin*(x-x0))^2+gamma^2)', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.05  1e2  spin  fPara(1)  hPara(1)]; % [wc g spin x0 gamma]
opts.Lower = [0    0   -Inf  fPara(2)  hPara(2)];
opts.Upper = [1    1    Inf  fPara(3)  hPara(3)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Resonant frequency fit' );
h = plot( fitresult, xData, yData );
legend( h, 'f0 vs. H0', 'Weak coupling fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'H0', 'Interpreter', 'none' );
ylabel( 'f0', 'Interpreter', 'none' );
grid on


