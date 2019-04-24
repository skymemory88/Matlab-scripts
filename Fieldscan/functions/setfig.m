function hfig = setfig(nfig)
% Function to setup the figure size and position
% INPUT:
% nfig      figure number
% OUTPUT:
% hfig      figure handle

hfig = figure(nfig);
clf
pos = get(hfig,'position');
set(hfig,'position',[pos(1:2) 600 500])


end