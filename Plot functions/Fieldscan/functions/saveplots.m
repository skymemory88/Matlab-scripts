function saveplots(hfig,figdir,figname)
% Function to save the figures to a specified directory and in different
% file formats
% INPUT:
% hfig      figure handle
% figdir    figure directory
% figname   file name to save figure to

    curdir = cd;
    cd(figdir)

    saveas(figure(hfig),[figname '.fig'],'fig');
%     print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
%     print(figure(hfig),[figname  '.png'],'-dpng','-r600');
    print2eps(figname,hfig)
    [~,~] = eps2xxx([figname '.eps'],{'jpeg','pdf'});

    disp(['Figure ' figname ' saved to '])
    disp(cd)
    cd(curdir)
end