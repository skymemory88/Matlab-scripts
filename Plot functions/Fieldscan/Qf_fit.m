function [fitresult, gof] = Qf_fit(H0, Q0, init, bdl, bdh, plt)
%CREATEFIT(H0,F0)
%  Create a fit.
%
%  Data for 'wk_kap_fit' fit:
%      X Input : H0
%      Y Output: Q0
%      hPara: [init_h, hBound_l, hBound_h]
%      fPara: [init_k, fBound_l, fBound_h]
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'wk_kap_fit'.
[xData, yData] = prepareCurveData( H0, Q0 );

% Set up fittype and options.
ft = fittype( @(gamma, gc, kappa, slope, x0, wc, x)...
    ((slope*(x-x0)).^2 + gamma^2)*wc ./ (2*gc^2*gamma + kappa*( (slope*(x-x0)).^2 + gamma^2)),...
    'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'gamma', 'gc', 'kappa', 'slope', 'x0', 'wc'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [init(1) init(2)  init(3)  init(4)  init(5)  init(6)]; % [slope x0 gamma gc kappa]
opts.Lower = [bdl(1)  bdl(2)  bdl(3)  bdl(4)  bdl(5)  bdl(6)];
opts.Upper = [bdh(1)  bdh(2)  bdh(3)  bdh(4)  bdh(5)  bdh(6)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
if plt
    figure( 'Name', 'Quality Factor fit' );
    h = plot( fitresult, xData, yData );
    legend( h, 'H0 vs. Qf', 'Q-fit', 'Location', 'SouthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'Magnetic field (T)', 'Interpreter', 'none' );
    ylabel( 'Quality Factor', 'Interpreter', 'none' );
    grid off
end


