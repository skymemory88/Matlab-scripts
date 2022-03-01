function [fitresult, gof] = iptopt(varargin)
%CREATEFIT(omega,S11)
%  S11 input-output fit:
%      X Input : frequency
%      Y Output: |S11|
%      coef : fixed parameter
%      param: fitting parameters (chi_off, w0 , Gc, gamma, attenuation) with upper and lower bound
%      aux_param: fixed parameters of an auxilliary field (for fit_mode 2 and 3)
%      weight: weight distribution among data pointss
%      fit_mode: fit model choice (0-3)
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
num_res = 6; % number of total spin resonances

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
opts.StartPoint = [P0(1) P0(2) P0(3) P0(4) P0(5) P0(6)  P0(7)]; % {'omega', 'Gc', 'gma', 'attn', 'xr'/'kpe', 'xi'/'kpi', 'wc'/'Gc1' }
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
    case 1 % fit with fixed cavity parameters
%         opts.Lower(5) = -Inf; % ['xi']
%         opts.Upper(5) =  Inf;
        opts.StartPoint(6:7) = []; % ['Gc1', 'Gc2']
        opts.Lower(6:7) = [];
        opts.Upper(6:7) = [];
        
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, xr, xi, kpe, kpi, wc, attn, x)...
            ( abs(1 + 2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma)) ).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients', {'omega', 'Gc', 'gma', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc','attn'});        
        % Fit model to data.
        [fitresult, gof] = fit(xData, yData, ft, opts, 'problem', {kpe, kpi, wc, attn});
    case 2 % fit with fixed cavity paramters and a single auxilliary spin resonance (auxilliary coupling strength loose)
%         nbr = mod(idx-1-1,num_res)+1;
        nbr = mod(idx+1-1,num_res)+1;
        w0 = aux_param.f0(nbr);
        
        opts.StartPoint(7) = aux_param.gc(nbr); % ['Gc2']
        opts.StartPoint(5:6) = []; % ['xi', 'Gc1']
        opts.Lower(5:6) = [];
        opts.Upper(5:6) = []; 
        
        % Set up fittype and options.
        ft = fittype( @(omega, Gc, gma, xr, Gc1, kpe, kpi, wc, attn, w0, x) (abs(1+2*kpe./...
            (1i*(x-wc+xr) - (kpe+kpi) + Gc^2./(1i*(x-omega)-gma) + Gc2^2./(1i*(x-w0)-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'},'coefficients', {'omega','Gc','gma','xr','Gc2'},...
            'problem', {'kpe','kpi','wc','w0','attn'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem', {kpe, kpi, wc, attn, w0});
    case 3 % fit with fixed cavity paramters and two auxilliary spin resonance (both auxilliary coupling strengths loose)
        lft = mod(idx-1-1, num_res)+1;
        rgt = mod(idx+1-1, num_res)+1;
        w1 = aux_param.f0(lft);
        w2 = aux_param.f0(rgt);

        opts.StartPoint(6) = aux_param.gc(lft); % ['Gc1']
        opts.StartPoint(7) = aux_param.gc(rgt); % ['Gc2']

        opts.StartPoint(5) = []; % ['xi']
        opts.Lower(5) = [];
        opts.Upper(5) = [];
        
        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, Gc1, Gc2, kpe, kpi, wc, attn, w1, w2, x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi) + Gc^2./(1i*(x-omega)-gma)...
            + Gc1^2./(1i*(x-w1)-gma) + Gc2^2./(1i*(x-w2)-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'Gc1', 'Gc2'},...
            'problem', {'kpe','kpi','wc','w1','w2','attn'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, attn, w1, w2});
    case 4 % fit with fixed cavity paramters and two auxilliary spin resonance (both auxilliary coupling strengths fixed)
        lft = mod(idx-1-1, num_res)+1;
        rgt = mod(idx+1-1, num_res)+1;
        w1 = aux_param.f0(lft);
        auxGc1 = aux_param.gc(lft);
        w2 = aux_param.f0(rgt);
        auxGc2 = aux_param.gc(rgt);

        opts.StartPoint(5) = 0; % ['xi']
        opts.Lower(5) = -Inf;
        opts.Upper(5) =  Inf;

        opts.StartPoint(6:7) = []; % ['Gc1', 'Gc2']
        opts.Lower(6:7) = [];
        opts.Upper(6:7) = [];
        
        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, kpe, kpi, wc, attn, w1, w2, auxGc1, auxGc2,x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma)...
            + auxGc1^2./(1i*(x-w1)-gma) + auxGc2^2./(1i*(x-w2)-gma))).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega','Gc','gma','xr','xi'},...
            'problem', {'kpe','kpi','wc','attn','w1','w2','auxGc1','auxGc2'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, attn, w1, w2, auxGc1, auxGc2});
    case 5 % fit with fixed cavity paramters and all five auxilliary spin resonances (all auxilliary coupling strengths fixed)
        aux_spin = 0;
        sr = 1:num_res;
        sr = sr(sr ~= idx); % remove the fitting target spin resonance
        for ii = sr
            w = aux_param.f0(ii);
            auxGc = aux_param.gc(ii);
            aux_spin = aux_spin + auxGc^2./(1i*(x-w)-2e-2);
        end

        opts.StartPoint(5) = 0; % ['xi']
        opts.Lower(5) = -Inf;
        opts.Upper(5) =  Inf;
        
        opts.StartPoint(6:7) = []; % ['wc','Gc1']
        opts.Lower(6:7) = [];
        opts.Upper(6:7) = [];

        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, kpe, kpi, wc, attn, aux_spin, x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma) + aux_spin )).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'xi'},...
            'problem', {'kpe','kpi','wc','attn','aux_spin'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts, 'problem',{kpe, kpi, wc, attn, aux_spin});
    case 6 % fit with fixed cavity paramters and all five auxilliary spin resonances (N.N. auxilliary coupling strengths loose)
        aux_spin = 0;
        lft = mod(idx-1-1,num_res)+1;
        rgt = mod(idx+1-1,num_res)+1;
        sr = 1:num_res;
        sr([lft idx rgt]) = []; % remove the fitting target spin resonance and it's nearest neighbours
        for ii = sr
            w = aux_param.f0(ii);
            auxGc = aux_param.gc(ii);
            aux_spin = aux_spin + auxGc^2./(1i*(x-w)-2e-2);
        end
        w1 = aux_param.f0(lft);
        w2 = aux_param.f0(rgt);

        % Set up fittype and options
        ft = fittype( @(omega, Gc, gma, xr, xi, Gc1, Gc2, kpe, kpi, wc, attn, w1, w2, aux_spin, x)...
            (abs(1+2*kpe./(1i*(x-wc+xr) - (kpe+kpi+xi) + Gc^2./(1i*(x-omega)-gma)...
            + Gc1^2./(1i*(x-w1)-gma) + Gc2^2./(1i*(x-w2)-gma) + aux_spin)).*attn),...
            'independent', {'x'}, 'dependent', {'y'}, 'coefficients',{'omega', 'Gc', 'gma', 'xr', 'xi', 'Gc1', 'Gc2'},...
            'problem', {'kpe','kpi','wc','attn','w1','w2','aux_spin'});
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts,...
            'problem',{kpe, kpi, wc, attn, w1, w2, aux_spin});
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
