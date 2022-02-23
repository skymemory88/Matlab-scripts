function [fitresult, gof] = iptopt(x, y, coef, param, aux_param, weight, fit_mode)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      coef : fixed parameter
%      P0: Initial guesses for fitting parameters (chi_off, w0 , Gc, gamma, attenuation)
%      Y Output: |S11|
%      Plt (true of false): Plot the fit along with data
%  Output:
%      fitresult : a fit object representing the S11 parameter.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
% Fit: 'S11_iptopt_fit'.
kpe = coef(1); % external dissipation rate
kpi = coef(2); % internal dissipation rate
wc = coef(3); % cavity resonance frequency
field = coef(4); % current magnetic field
plt = coef(5); % plotting option

P0 = param(1,:); % starting values of the fitting parameters
low = param(2,:); % lower limit of the fitting parameters
upr = param(3,:); % upper limit of the fitting parameters

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)  P0(7)]; % {'omega', 'Gc', 'gma', 'attn', 'xr'/'kpe', 'xi'/'kpi', 'wc'/'Gc1' }
opts.Lower = [low(1) low(2) low(3) low(4) low(5) low(6) low(7)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5) upr(6) upr(7)];
opts.Weights = weight;

[xData, yData] = prepareCurveData( x, y );
switch fit_mode
    case 0 % fit all parameters with single spin resonance
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, attn, kpe, kpi, wc, x) (abs(1 + 2*kpe./...
            (1i*(x-wc) - (kpi + kpe) + Gc^2./(1i*(x-omega )-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'attn', 'kpe', 'kpi', 'wc'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts);
    case 1 % fit with fixed cavity parameters
        opts.StartPoint(7) = []; % {'xr', 'xi'}
        opts.Lower(7) = [];
        opts.Upper(7) = [];
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, attn, xr, xi, kpe, kpi, wc, x)...
            ( abs(1 + 2*kpe./(1i*(x - wc + xr) - (kpi + kpe + xi) + Gc^2./(1i*(x-omega)-gma)) ).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients', {'omega', 'Gc', 'gma', 'attn', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc'});        
        % Fit model to data.
        [fitresult, gof] = fit(xData, yData, ft, opts, 'problem', {kpe, kpi, wc});
    case 2 % fit with fixed cavity paramters and a single auxilliary spin resonance (loose auxilliary coupling strength)
        w2 = aux_param.f0(2);
        opts.StartPoint(5:6) = 0; % {'xr', 'xi'}
        opts.Lower(5:6) = 0;
        opts.Upper(5:6) = 0;
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, attn, xr, xi, Gc1, kpe, kpi, wc, w2, x) (abs(1+2*kpe./...
            (1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma) + Gc1^2./(1i*(x-w2)-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'},'coefficients', {'omega', 'Gc', 'gma','attn', 'xr', 'xi', 'Gc1'},...
            'problem', {'kpe','kpi','wc','w2'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {kpe, kpi, wc, w2});
    case 3 % fit with fixed cavity paramters and two auxilliary spin resonance (fixed auxilliary coupling strengths)
        axGc1 = aux_param.gc(1);
        w1 = aux_param.f0(1);
        axGc2 = aux_param.gc(2);
        w2 = aux_param.f0(2);
        
        opts.StartPoint(5:6) = 0; % {'xr', 'xi'}
        opts.Lower(5:6) = 0;
        opts.Upper(5:6) = 0;       
        opts.StartPoint(7) = []; % {'xr', 'xi'}
        opts.Lower(7) = [];
        opts.Upper(7) = [];
        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, attn, xr, xi, kpe, kpi, wc, axGc1, w1, axGc2, w2, x)...
            (abs(1 + 2*kpe./(1i*(x - wc + xr) - (kpi + kpe + xi) + Gc^2./(1i*(x-omega)-gma)...
            + axGc1.^2./(1i*(x-w1)-gma) + axGc2.^2./(1i*(x-w2)-gma) )).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'attn', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc','axGc1','w1','axGc2','w2'});        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, axGc1, w1, axGc2, w2});
end

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
