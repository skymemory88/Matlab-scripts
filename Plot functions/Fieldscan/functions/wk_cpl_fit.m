function [fitresult, gof] = wk_cpl_fit(H0, f0, init, bdl, bdh, weight, plt)
%CREATEFIT(H0,F0)
%  Strong coupling fit of resonant frequency.
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
ft = fittype( 'wc - g^2*slope*(x-x0)/((slope*(x-x0))^2 + gamma^2)', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [init(1) init(2)  init(3)  init(4)  init(5)]; % [g  gamma slope wc x0]
opts.Lower = [bdl(1)  bdl(2)  bdl(3)  bdl(4)  bdl(5)];
opts.Upper = [bdh(1)  bdh(2)  bdh(3)  bdh(4)  bdh(5)];
opts.Weights = weight;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts);

% Plot fit with data
if plt
    figure( 'Name', 'Resonant frequency fit' );
    h = plot( fitresult, xData, yData );
    legend( h, 'f0 vs. H0', 'Weak coupling fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'H0', 'Interpreter', 'none' );
    ylabel( 'f0', 'Interpreter', 'none' );
    grid off
end


