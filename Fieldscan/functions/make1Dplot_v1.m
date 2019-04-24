function varargout = make1Dplot_v1(hfig,out,varx,vary,plotopt)
% Function to plot the data in the out structure
%
% INPUT:
% hfig:             handle of figure to which to plot data
% varx:             name of field to use as the x-axis
% vary:             name of field to use as the y-axis
% plotopt:          additional parameters for a nice plot
% plotopt.trans:    [0] or 1 define a mapping of parameters x,y,dx,dy, of
%                   the form:
%   plotopt.transx:   x' = x;
%   plotopt.transdx:  dx' = dx;
%   plotopt.transy:   y' = 1./y;
%   plotopt.transdy:  dy' = dy./y.^2
% plotopt.axes:     [0] or 1 to plot into axes defined by handle plotopt.axes
% plotopt.labels:   'x','y' or 'xy' - define which axes to label
% plotopt.ftsz:     Define fontsize of the labels
%
% OUTPUT:
% h:                handles in a cell structure to the plotted data
% ax:               handle to the plot axes

%% Define figure position on screen and axes

figure(hfig)
set(hfig,   'position',[1 1 600 500],...
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

% Determine if fontsize are defined or use defaults
if ~isfield(plotopt,'ftsz')
    plotopt.ftsz = 12;
end

box on
hold on

%% Plot data to figure

for n = 1:length(out)
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
            eval(['x = ' plotopt.transx ';'])
            eval(['dx = ' plotopt.transdx ';'])
            eval(['y = ' plotopt.transy ';'])
            eval(['dy = ' plotopt.transdy ';'])
        end
    end
           
   
    h{n} = ploterr(x,y,dx,dy,'hhxy',0);
    set(h{n},       'Color',plotopt.col,...
                    'linewidth',plotopt.lnwd)
    set(h{n}(1),    'MarkerSize',plotopt.mksz,...
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
    xlabel(varx,'interpreter','latex','fontsize',plotopt.ftsz)
    ylabel(vary,'interpreter','latex','fontsize',plotopt.ftsz)
end

set(ax,'fontsize',plotopt.ftsz)

hold off

% Define output
if nargout == 1
    varargout = h;
end
if nargout == 1
    varargout = {h, ax};
end


end
