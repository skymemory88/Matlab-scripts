function [fitresult, gof] = gcFit(x, y, temp)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% meV2J = 1.602217e-22; % [J/meV]
kB = 0.08617; % [meV/K]
% hbar = 1.055E-34/meV2J*1e9; % Reduced Planck constant [meV/GHz]

% Set up fittype and options.
% ft = fittype( 'a*x*(exp(-hbar*b*x/kB/temp)-exp(-hbar*c*x/kB/temp))', 'independent', {'x'}, 'dependent', {'y'} , ...
%     'coefficient',{'a','b','c'}, 'problem', {'temp', 'kB', 'hbar'});

ft = fittype( 'a*x*(exp(-c*x/temp)-exp(-d*x/temp)/b)', 'independent', {'x'}, 'dependent', {'y'} , ...
    'coefficient',{'a','b','c','d'},'problem',{'temp'});
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.1 0.1 0.1 0.1];
opts.Lower(1) = -Inf;
opts.Upper(1) = Inf;
opts.Lower(2:4) = 0;
opts.Upper(2:4) = Inf;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts,'problem',{temp});
% [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{temp, kB, hbar});

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid off


