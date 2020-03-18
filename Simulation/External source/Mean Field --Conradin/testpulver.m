for n=1:100
    [fpar,sp,Nq]=chiqw_pulver(1.5,0.01,1+n*0.25);
    figure
    plot(sp)
    spt(n)=sp(1);
    for m=2:Nq
        spt(n)=spt(n)+sp(m);
    end
end
%save pulverwerte spt
temp=[1.25:0.25:26];
mapplot(spt,temp);
shading interp

load pulverwerte spt
for n=1:100
    [spf,fp]=fits(spt(n),'lorz',[1000,0.9,0.01,0],[1 1 1 0]);
    plot(spf)
    breite(n)=fp.pvals(3);
    pos(n)=fp.pvals(2);
end
spb=spec1d(temp,abs(breite),0.001*abs(breite));
plot(spb)
figure
spp=spec1d(temp,abs(pos),0.0000001*abs(pos));
plot(spp)