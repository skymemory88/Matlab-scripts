function sendtomfit(s1,varargin)
%
% function sendtomfit
%
% SPEC1D/PLOT Send spec1d object to mfit
%
% PB 23.07.12

if nargin == 1
    selected = ones(size(getfield(s1,'x')));
    p = [];
    fixed = [];
end
if nargin == 2
    selected = varargin{1};
    p = [];
    fixed = [];
end
if nargin == 3
    selected = varargin{1};
    p = varargin{2};
    fixed = [];
end
if nargin == 4
    selected = varargin{1};
    p = varargin{2};
    fixed = varargin{3};
end


[x y e] = extract(s1);


tomfit(x,y,e,selected,p,fixed);

fig = findobj('tag','mf_DataWindow');

figure(fig);
set(fig,'color','w')
box on

xlim([0.95*min(x) 1.05*max(x)])
ylim([0.95*min(y) 1.05*max(y)])

set(gca,'xcolor','k')
set(gca,'ycolor','k')

try
    xlabel(getfield(s1,'x_label'),'color','k')
    ylabel(getfield(s1,'y_label'),'color','k')
catch ME
end



end