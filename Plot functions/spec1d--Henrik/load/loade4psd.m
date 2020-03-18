function [d,sd,s]=loade4psd(path,files,winx,winy,xstr,nmstr,flag)
[d,sd]=e4data([path,'\e4',num2str(files),'.asc']);
dd = hdfread([path,'\0',num2str(files),'.nxs'], '/entry1/data1/DATA', 'Index', {[1  1  1],[1  1  1],[size(d,1)  256  256]});
dd=double(dd);
if nargin==7
%    winx=1:256;
%    winy=1:256;
    pcolor(squeeze(sum(dd(:,winy,winx))));grid off;shading interp
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
