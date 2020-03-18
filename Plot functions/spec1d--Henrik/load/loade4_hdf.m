function [d,sd,s]=loade4psd(path,files,winx,winy,xstr,nmstr,flag)
[d,sd]=e4data([path,'\e4',num2str(files),'.asc']);
dd = hdfread([path,'\',num2str(files),'_LAMP.hdf'], '/entry1/data1/DATA', 'Index', {[1  1  1],[1  1  1],[size(d,1)  256  256]});
dd=double(dd);
if nargin==7
    switch flag
        case 'all'
            winx=1:256;
            winy=1:256;
            pcolor(squeeze(sum(dd(:,winy,winx))));grid off;shading interp
        case 'win'
            pcolor(squeeze(sum(dd(:,winy,winx))));grid off;shading interp
        case 'movie'
            for n=1:size(dd,1);
                dimg=log(squeeze(dd(n,winy,winx)));
                M(n).cdata=uint8(round(dimg/max(max(dimg))*256));
                M(n).colormap=hot(256);
            end
            implay(M)
        case 'montage'
            m=ceil(sqrt(size(dd,1)));
            for n=1:size(dd,1)
                subplot(m,m,n,'align')
                dimg=(squeeze(dd(n,winy,winx)));
                h=pcolor(dimg);grid off;shading interp
                set(gca,'xticklabel','')    
                set(gca,'yticklabel','')
                text(2,15,num2str(n),'color','w','fontsize',10)
            end
        case 's1d'
    end
end

ds=squeeze(sum(sum(dd(:,winy,winx),3),2));
no=strmatch(xstr,sd);
if nmstr == 'none'
    nm=1;
else
    nm=strmatch(nmstr,sd);
end
s=spec1d(d(:,no),ds,sqrt(ds))*nm;

return
