function varargout = makeplot_v1(hfig,out,varx,vary,plotoptin)
% Function to plot the data in the out structure
%
% INPUT:
% hfig:                     handle of figure to which to plot data
% varx:                     name of field to use as the x-axis
% vary:                     name of field to use as the y-axis
% plotopt:                  additional parameters for a nice plot
% plotopt.trans:            [0] or 1 define a mapping of parameters x,y,dx,dy, of
%                           the form:
%                           plotopt.transx:   x' = x;
%                           plotopt.transdx:  dx' = dx;
%                           plotopt.transy:   y' = 1./y;
%                           plotopt.transdy:  dy' = dy./y.^2
% plotopt.axes:             [0] or 1 to plot into axes defined by handle plotopt.axes
% plotopt.labels:           'x','y' or 'xy' - define which axes to label
% plotopt.ftsz:             Define fontsize of the labels
% plotopt.col:              Define colour of points, this can either be a [r g b]
%                           array or input for "setcolours.m" function, ie 'jet'
% plotopt.xlabel/ylabel:    set labels for axes
% plotopt.title:            set title for axes
% plotopt.noerror:          [0] or 1 Remove error bars from the plot
%                   
%
% OUTPUT:
% h:                handles in a cell structure to the plotted data
% ax:               handle to the plot axes

%% Define default parameters for the plot =================================
plotopt.col  = [0.2 0.2 0.7];                       % Plot colour
plotopt.lnwd = 1;                                   % Linewidth
plotopt.ftsz = 12;                                  % Text fontsize
plotopt.mksz = 5;                                   % Marker size
plotopt.noerror = 0;                                % Plot error bars if 0

if nargin > 4
    % Overwrite default settings
    fnames = fieldnames(plotoptin);
    for m = 1:length(fnames)
        plotopt.(fnames{m}) = plotoptin.(fnames{m});
    end
end


%% Define figure position on screen and axes

figure(hfig)
try
    pos = get(hfig,'position');
catch ME
    pos = [1 1 0 0];
end
set(hfig,   'position',[pos(1:2) 600 500],...
            'tag','fig_ac_suscept')

if isfield(plotopt,'axes')
    % Plot to existing axes defined by user    
    ax = plotopt.axes;
else
    % Determine whether or not axes already exist in the figure window
    ax = findobj(hfig,'tag','ax_ac_suscept');
    
    if isempty(ax)
        ax = axes('position',[0.15 0.15 0.75 0.75],'tag','ax_ac_suscept');
    end
end
axes(ax);


% Set default colour scheme
if isfield(plotopt,'col')
    if ischar(plotopt.col)
        plotcol_flag = plotopt.col;
    end
else
    plotopt.col = [0.2 0.2 0.7];
end


box on
hold on

%% Plot data to figure

for n = 1:length(out)
    
    if exist('plotcol_flag','var')
        if length(out) < 5
            plotopt.col = setcolours([],num2str(n));
        else
            plotopt.col = setcolours(n/length(out),plotcol_flag);
        end
    end
        
    % Determine if errors have been calculated
    if ~isfield(out(n).data,['d' varx])
        out(n).data.(['d' varx]) = out(n).data.(varx).*0;
    end
    if ~isfield(out(n).data,['d' vary])
        out(n).data.(['d' vary]) = out(n).data.(vary).*0;
    end

    x = out(n).data.(varx);
    y = out(n).data.(vary);
    dx = out(n).data.(['d' varx]);
    dy = out(n).data.(['d' vary]);
    
    % Transform x,y,dx,dy if necessary
    if isfield(plotopt,'trans')
        if plotopt.trans == 1
            if isfield(plotopt,'transx')
                eval(['x = ' plotopt.transx ';'])
            end
            if isfield(plotopt,'transdx')
                eval(['dx = ' plotopt.transdx ';'])
            end
            if isfield(plotopt,'transy')
                eval(['y = ' plotopt.transy ';'])
            end
            if isfield(plotopt,'transdy')
                eval(['dy = ' plotopt.transdy ';'])
            end
        end
    end

    if plotopt.noerror == 0
        h(:,n) = ploterr(x,y,dx,dy,'hhxy',0);
    else
        h(:,n) = ploterr(x,y,[],[],'hhxy',0);
    end
    if plotopt.lnwd > 0
        set(h(:,n), 'Color',plotopt.col,...
                    'linewidth',plotopt.lnwd)
    else
        set(h(:,n), 'Color',plotopt.col)
    end
    set(h(1,n),    'MarkerSize',plotopt.mksz,...
                    'MarkerEdgeColor',plotopt.col,...
                    'MarkerFaceColor',plotopt.col,...
                    'Marker','o',...
                    'LineStyle','none')
end

% Label plot axes if necessary
if isfield(plotopt,'labels')
    if strcmpi(plotopt.labels,'x') || strcmpi(plotopt.labels,'xy')        
        xlabel(varx,'interpreter','latex','fontsize',plotopt.ftsz)
    end
    if strcmpi(plotopt.labels,'y') || strcmpi(plotopt.labels,'xy')        
        ylabel(vary,'interpreter','latex','fontsize',plotopt.ftsz)
    end
else
    if isfield(plotopt,'xlabel'), varx = plotopt.xlabel; end
    if isfield(plotopt,'ylabel'), vary = plotopt.ylabel; end
    xlabel(varx,'interpreter','latex','fontsize',plotopt.ftsz)
    ylabel(vary,'interpreter','latex','fontsize',plotopt.ftsz)
end

% Add title to plot
if isfield(plotopt,'title')
    title(plotopt.title,'interpreter','latex','fontsize',plotopt.ftsz)
end

set(ax,'fontsize',plotopt.ftsz)

hold off

% Define output
if nargout == 1
    varargout{1} = h;
end
if nargout == 2
    varargout{1} = h;
    varargout{2} = ax;
end


end
