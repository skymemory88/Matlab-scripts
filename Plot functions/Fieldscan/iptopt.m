function [fitresult, gof] = iptopt(varargin)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      Y Output: |S11|
%      coef : fixed parameter
%      param: fitting parameters (chi_off, w0 , Gc, gamma, attenuation) with upper and lower bound
%      aux_param: fixed parameters of an auxilliary field (for fit_mode 2 and 3)
%      weight: weight distribution among data pointss
%      ioFit: fit function setting
%         init: default to fit mode 0 (for empty cavity fitting)
%         mode: fit model choice (0-6)
%         tarIdx: spin resonance frequency index (1-6)
%  Output:
%      fitresult : a fit object representing the S11 parameter.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
% Fit: 'S11_iptopt_fit'.
x = varargin{1};
y = varargin{2};
coef = varargin{3};
param = varargin{4};
aux_param = varargin{5};
weight = varargin{6};
ioFit = varargin{7};

if ioFit.init == true; ioFit.mode = 0; end

idx = ioFit.tarIdx;
num_res = ioFit.total-1; % total number of total spin resonances
if  idx > 2 && idx < num_res
    NN1 = mod(idx+1-1, num_res)+1; % Nearest neighbour (right)
    NN2 = mod(idx-1-1, num_res)+1; % Nearest neighbour (left)
elseif idx == num_res % right most spin resonance on PM side
    NN1 = mod(idx-1-1, num_res)+1; % Nearest neighbour
    NN2 = mod(idx-2-1, num_res)+1; % Next Nearest neighbour
elseif idx == 2 % left most spin resonance on PM side
    NN1 = mod(idx+1-1, num_res)+1; % Nearest neighbour
    NN2 = mod(idx+2-1, num_res)+1; % Next Nearest neighbour
elseif idx == 1 % spin resonance on FM side
    % do nothing (no neighbours on FM side)
end

kpe = coef(1); % external dissipation rate
kpi = coef(2); % internal dissipation rate
wc = coef(3); % cavity resonance frequency
attn = coef(4); % attenuation rate
field = coef(5); % current magnetic field
plt = coef(6); % plotting option

P0 = param(1,:); % starting values of the fitting parameters
low = param(2,:); % lower limit of the fitting parameters
upr = param(3,:); % upper limit of the fitting parameters

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
% {'omega', 'Gc', 'gma', 'xr/kpe', 'xi/kpi', 'wc/Gc1', 'attn/Gc2'}
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)  P0(7)];
opts.Lower = [low(1) low(2) low(3) low(4) low(5) low(6) low(7)];
opts.Upper = [upr(1) upr(2) upr(3) upr(4) upr(5) upr(6) upr(7)];
opts.Weights = weight;

[xData, yData] = prepareCurveData( x, y );
switch ioFit.mode
    case 0 % fit all parameters with single spin resonance
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, kpe, kpi, wc, attn, x) (abs(1 + 2*kpe./...
            (1i*(x-wc) - (kpi+kpe) + Gc^2./(1i*(x-omega)-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'kpe', 'kpi', 'wc', 'attn'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts);
    case 1
        % fit with fixed cavity parameters
        opts.StartPoint(6:7) = []; % ['Gc1', 'Gc2']
        opts.Lower(6:7) = [];
        opts.Upper(6:7) = [];
        
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, xr, xi, kpe, kpi, wc, attn, x)...
            ( abs(1 + 2*kpe./(1i*(x-wc) - (kpe+kpi) + 1i*(xr+1i*xi) + Gc^2./(1i*(x-omega)-gma)) ).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients', {'omega', 'Gc', 'gma', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc','attn'});        
        % Fit model to data.
        [fitresult, gof] = fit(xData, yData, ft, opts, 'problem', {kpe, kpi, wc, attn});
    case 2 
        % fit to one neighbouring spin resonances with fixed cavity paramters
        w1 = aux_param.f0(NN1);
        gma1 = aux_param.gma(NN1);
        opts.StartPoint(6) = aux_param.gc(NN1); % ['Gc1']
        
        opts.StartPoint(7) = w1; % ['w1']
%         opts.Lower(7) = w1;
%         opts.Upper(7) = w1; 
        opts.Lower(7) = -Inf;
        opts.Upper(7) = Inf; 

%         opts.StartPoint(7) = []; % ['Gc2']
%         opts.Lower(7) = [];
%         opts.Upper(7) = []; 
        
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, xr, xi, Gc1, w1, kpe, kpi, wc, attn, gma1, x) (abs(1+2*kpe./...
            (1i*(x-wc) - (kpe+kpi) + 1i*(xr+1i*xi) + Gc^2./(1i*(x-omega)-gma) + Gc1^2./(1i*(x-w1)-gma1))).*attn),...
            'independent', {'x'}, 'dependent', {'y'},'coefficients', {'omega','Gc','gma','xr','xi','Gc1','w1'},...
            'problem', {'kpe','kpi','wc','attn','gma1'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {kpe, kpi, wc, attn, gma1});
    case 3
        % fit with fixed cavity paramters and two auxilliary spin resonance (auxilliary coupling strengths fixed, spin linewidth loose)
        w1 = aux_param.f0(NN1);
        aGc1 = aux_param.gc(NN1);
        gma1 = aux_param.gma(NN1);
        
        w2 = aux_param.f0(NN2);
        aGc2 = aux_param.gc(NN2);
        gma2 = aux_param.gma(NN2);

        opts.StartPoint(6) = gma1; % ['gma1']
        opts.StartPoint(7) = gma2; % ['gma2']
%         opts.StartPoint(6:7) = []; % ['Gc1', 'Gc2']
%         opts.Lower(6:7) = [];
%         opts.Upper(6:7) = [];
        
        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, gma1, gma2, kpe, kpi, wc, attn, w1, w2, aGc1, aGc2, x)...
            (abs(1+2*kpe./(1i*(x-wc) - (kpe+kpi) + 1i*(xr+1i*xi) + Gc^2./(1i*(x-omega)-gma)...
            + aGc1^2./(1i*(x-w1)-gma1) + aGc2^2./(1i*(x-w2)-gma2)) ).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega','Gc','gma','xr','xi','gma1','gma2'},...
            'problem', {'kpe','kpi','wc','attn','w1','w2','aGc1','aGc2'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, attn, w1, w2, aGc1, aGc2});
    case 4 
        % fit with fixed cavity paramters and two auxilliary spin resonance (auxilliary coupling strengths loose, spin linewidths fixed)
        w1 = aux_param.f0(NN1);
        gma1 = aux_param.gma(NN1);
        w2 = aux_param.f0(NN2);
        gma2 = aux_param.gma(NN2);

        opts.StartPoint(6) = aux_param.gc(NN1); % ['Gc1']
        opts.StartPoint(7) = aux_param.gc(NN2); % ['Gc2']
%         opts.StartPoint(6) = aux_param.gma(NN1); % ['gma1']
%         opts.StartPoint(7) = aux_param.gma(NN2); % ['gma2']
        
        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, Gc1, Gc2, kpe, kpi, wc, w1, w2, attn, gma1, gma2, x)...
            (abs(1+2*kpe./(1i*(x-wc) - (kpe+kpi) + 1i*(xr+1i*xi) + Gc^2./(1i*(x-omega)-gma)...
            + Gc1^2./(1i*(x-w1)-gma1) + Gc2^2./(1i*(x-w2)-gma2))).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'xi', 'Gc1', 'Gc2'},...
            'problem', {'kpe','kpi','wc','w1','w2','attn','gma1','gma2'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, w1, w2, attn, gma1, gma2});
    case 5 
        % fit with fixed cavity paramters and all auxilliary spin resonances (all auxilliary coupling strengths fixed)
        aux_spin = 0;
        sr = 1:num_res;
        sr = sr(sr ~= idx); % remove the fitting target spin resonance
        for ii = sr
            w = aux_param.f0(ii);
            aGc = aux_param.gc(ii);
            gma = aux_param.gma(ii);
            aux_spin = aux_spin + aGc^2./(1i*(x-w)-gma);
        end
        
        opts.StartPoint(6:7) = []; % ['Gc1','Gc2']
        opts.Lower(6:7) = [];
        opts.Upper(6:7) = [];

        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, kpe, kpi, wc, attn, aux_spin, x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma) + aux_spin )).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc','attn','aux_spin'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, attn, aux_spin});
    case 6 
        % fit with fixed cavity paramters and all five auxilliary spin resonances (N.N. auxilliary coupling strengths loose)
        aux_spin = 0;
        sr = 1:num_res;
        sr([NN1 idx NN2]) = []; % remove the fitting target spin resonance and it's nearest neighbours
        for ii = sr
            w = aux_param.f0(ii);
            aGc = aux_param.gc(ii);
            gma = aux_param.gma(ii);
            aux_spin = aux_spin + aGc^2./(1i*(x-w)-gma);
        end
        w1 = aux_param.f0(NN1);
        gma1 = aux_param.gma(NN1);
        w2 = aux_param.f0(NN2);
        gma2 = aux_param.gma(NN2);
        opts.StartPoint(6) = aux_param.gc(NN1); % ['Gc1']
        opts.StartPoint(7) = aux_param.gc(NN2); % ['Gc2']
%         opts.StartPoint(6) = opts.StartPoint(2); % ['Gc1']
%         opts.StartPoint(7) = opts.StartPoint(2); % ['Gc2']

        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, Gc1, Gc2, kpe, kpi, wc, attn, w1, w2, gma1, gma2, aux_spin, x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma)...
            + Gc1^2./(1i*(x-w1)-gma1) + Gc2^2./(1i*(x-w2)-gma2) + aux_spin)).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'xi', 'Gc1', 'Gc2'},...
            'problem', {'kpe','kpi','wc','attn','w1','w2','gma1','gma2','aux_spin'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {kpe, kpi, wc, attn, w1, w2, gma1, gma2, aux_spin});
end

if plt
    % Plot fit with data.
    figure( 'Name', sprintf('S11_inptopt_fit at B = %.3f T',field) );
    h = plot( fitresult, xData, yData );
    legend( h, 'Data', 'S11_fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    set(h, 'LineWidth', 1.5);
    % Label axes
    xlabel( 'Frequency (GHz)', 'Interpreter', 'none' );
    ylabel( '|S11|', 'Interpreter', 'none' );
    grid on
end
