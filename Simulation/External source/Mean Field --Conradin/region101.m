load werte/pgn pt ph

epsilon=0.002;
omega=[0.1:0.001:0.65];
H=[];
omgs=[];
dh=[-0.1:0.025:0];
dk=[-0.1:0.025:0];
for n=1:length(dh)
    for m=1:length(dk)
    S=chiqwn([1+dh(n) 0 1+dk(m)], omega, epsilon, ph(3)+0.001, pt(3));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega;
    S(nanind)=[];
    omegan(nanind)=[];
    s(n,m)=spec1d(omegan,S,sqrt(S));
    %[wmax,imax]=max(S);
    %sfold(n)=smooth(s(n,m),0.01);
    %plot(sfold)
    %[sf,f]=fits(sfold(n),'ngauss',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    %plot(sf)
    %pause(0.1)
    end
end
for nd=1:size(s,1)*size(s,2)
s(nd)=cut(s(nd),[0.1,1]);
end
%save werte/r101 s dh dk
for n=1:length(dh)
    for m=1:length(dk)
    sfold=smooth(s(n,m),0.01);
    y=getfield(sfold,'y');
    x=getfield(sfold,'x');
    ind=find(y>100)
    %plot(sfold)
    [sf,f]=fits(sfold,'ngauss',[100,x(ind(1)),0.05,100,x(ind(end)),0.05,0],[1 1 1 1 1 1 0]);
    e(n,m)=f.pvals(2);
    e2(n,m)=f.pvals(5);
    plot(sf)
    pause(0.1)
    end
end

mesh(e)