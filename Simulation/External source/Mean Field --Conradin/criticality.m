pt=[0:0.05:2.3];
ph=[];
h=[0:0.01:4.8];
for nt=1:length(pt)
    jz=3;
    jx=2;
    jy=0;
    jxyz=[];
    for nh=1:length(h)
        if jz>0.04
            jz=jz-0.01;
        end
        [jx,jy,jz]=remfn(h(nh),pt(nt),jx,jy,jz);
        jxyz=[jxyz;[jx,jy,jz]];
        h(nh);
    end
%     figure
%     plot(diff(jxyz(:,3)))
%     figure
%     plot(jxyz(:,3))
%     figure
%     plot(diff(jxyz(:,1)))
    [m,ind]=min(diff(jxyz(:,3)));
    ph=[ph,h(ind+1)];
end

plot(pt,ph,'.')

%save werte/pgn pt ph
load werte/pgn pt ph

epsilon=0.001;
omega=[0.025:0.0005:0.6];
omega2=[0.025:0.005:0.7];
H=[];
omgs=[];
for n=1:length(pt)
    S=chiqwn([1.95 0 0], omega, epsilon, ph(n), pt(n));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega;
    S(nanind)=[];
    omegan(nanind)=[];
    s(n)=spec1d(omegan,S,sqrt(S));
    [wmax,imax]=max(S);
    sfold(n)=smooth(s(n),0.01);
    %plot(sfold)
    [sf,f]=fits(sfold(n),'gaussskew',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    e(n)=f.pvals(2);
    plot(sf)
    pause(0.1)
    
    S=chiqwn([0.95 0 0], omega2, 10*epsilon, ph(n), pt(n));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega2;
    S(nanind)=[];
    omegan(nanind)=[];
    s2(n)=spec1d(omegan,S,sqrt(S));
    [wmax,imax]=max(S);
    sfold2(n)=smooth(s2(n),0.01);
    %plot(sfold)
    [sf,f]=fits(sfold2(n),'gaussskew',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    e2(n)=f.pvals(2);
    plot(sf)
    pause(0.1) 
%     %waitforbuttonpress()
%     pause(0.1)
%     
%     S=chiqwn([1 0 0], omega2, epsilon, ph(n)+0.01, pt(n));
%     nanind=find(isnan(S)|isinf(S));
%     omegan=omega2;
%     S(nanind)=[];
%     omegan(nanind)=[];
%     ind=find(S>100);
%     s2(n)=spec1d(omegan,S,sqrt(S));
%     [sf,f]=fits(cut(s2(n),[-0.02,0.02]+omegan(ind(end))),'lorz',[100,omegan(ind(end)),0.002,0],[1 1 1 0]);
%     e2(n)=f.pvals(2);
%     plot(sf)
%     %waitforbuttonpress()
%     pause(0.1)
end

figure
get(gcf)
set(gcf,'PaperPosition',[0.634517 6.34517 21/3 16])
set(gcf,'Position',[420 28 187*2 2*420])
axes('Position',[0.15,0.2/3,0.8,0.8/3])
h=plot(pt(1:35),e(1:35)./e2(1:35),'.')
ylabel('E_c/\Delta')
xlabel('T   K')
set(gca,'XLim',[0,2])
axes('Position',[0.15,(0.1+0.2)/3+0.8/3,0.8,0.8/3])
h2=plot(pt(1:35),e2(1:35),'.',pt(1:35),e(1:35),'.')
ylabel('E   meV')
set(gca,'XLim',[0,2])
axes('Position',[0.15,(2*0.1+0.2)/3+2*0.8/3,0.8,0.8/3])
h3=plot(pt,ph)
ylabel('H_c   Tesla')
set(gca,'XLim',[0,2])
legend(h2,['100';'200'])
print -depsc2 werte\gap_2.eps

%save werte/modes_2 e e2

%%%%% ohne hyp
pt=[0:0.05:2.3];
ph=[];
h=[0:0.01:4.8];
for nt=1:length(pt)
    jz=3;
    jx=2;
    jy=0;
    jxyz=[];
    for nh=1:length(h)
        if jz>0.04
            jz=jz-0.01;
        end
        [jx,jy,jz]=remfn(h(nh),pt(nt),jx,jy,jz);
        jxyz=[jxyz;[jx,jy,jz]];
        h(nh);
    end
%     figure
%     plot(diff(jxyz(:,3)))
%     figure
%     plot(jxyz(:,3))
%     figure
%     plot(diff(jxyz(:,1)))
    [m,ind]=min(diff(jxyz(:,3)));
    ph=[ph,h(ind+1)];
end

plot(pt,ph,'.')

%save werte/pgn pt ph
save werte/pgel pt ph

epsilon=0.001;
omega=[-0.3:0.005:0.6];
omega2=[-0.25:0.005:0.7];
H=[];
omgs=[];
for n=1:length(pt)
    S=chiqwn([1.95 0 0], omega, epsilon, ph(n), pt(n));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega;
    S(nanind)=[];
    omegan(nanind)=[];
    s(n)=spec1d(omegan,S,sqrt(S));
    [wmax,imax]=max(S);
    sfold(n)=smooth(s(n),0.05);
    %plot(sfold)
    [sf,f]=fits(sfold(n),'gaussskew',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    e(n)=f.pvals(2);
    plot(sf)
    pause(0.1)
    
    S=chiqwn([1 0 0], omega2, 10*epsilon, ph(n), pt(n));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega2;
    S(nanind)=[];
    omegan(nanind)=[];
    s2(n)=spec1d(omegan,S,sqrt(S));
    [wmax,imax]=max(S);
    sfold2(n)=smooth(s2(n),0.05);
    %plot(sfold)
    [sf,f]=fits(sfold2(n),'gaussskew',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    e2(n)=f.pvals(2);
    plot(sf)
    pause(0.1) 
%     %waitforbuttonpress()
%     pause(0.1)
%     
%     S=chiqwn([1 0 0], omega2, epsilon, ph(n)+0.01, pt(n));
%     nanind=find(isnan(S)|isinf(S));
%     omegan=omega2;
%     S(nanind)=[];
%     omegan(nanind)=[];
%     ind=find(S>100);
%     s2(n)=spec1d(omegan,S,sqrt(S));
%     [sf,f]=fits(cut(s2(n),[-0.02,0.02]+omegan(ind(end))),'lorz',[100,omegan(ind(end)),0.002,0],[1 1 1 0]);
%     e2(n)=f.pvals(2);
%     plot(sf)
%     %waitforbuttonpress()
%     pause(0.1)
end

figure
get(gcf)
set(gcf,'PaperPosition',[0.634517 6.34517 21/3 16])
set(gcf,'Position',[420 28 187*2 2*420])
axes('Position',[0.15,0.2/3,0.8,0.8/3])
h=plot(pt(1:35),e(1:35)./e2(1:35),'.')
ylabel('E_c/\Delta')
xlabel('T   K')
set(gca,'XLim',[0,2])
axes('Position',[0.15,(0.1+0.2)/3+0.8/3,0.8,0.8/3])
h2=plot(pt(1:35),e2(1:35),'.',pt(1:35),e(1:35),'.')
ylabel('E   meV')
set(gca,'XLim',[0,2])
axes('Position',[0.15,(2*0.1+0.2)/3+2*0.8/3,0.8,0.8/3])
h3=plot(pt,ph)
ylabel('H_c   Tesla')
set(gca,'XLim',[0,2])
legend(h2,['100';'101'])
print -depsc2 werte\gap_withouthyperfine.eps

%save werte/modes e e2
save werte/modesel e e2


h=[0:0.1:4.8];
dipolq([1.01 0 1;1.01 0 0;2.01 0 0;2 0 0],1,50)
for n=1:length(h)
    S=chiqwn([1.95 0 0], omega, epsilon, h(n), pt(2));
    nanind=find(isnan(S)|isinf(S));
    omegan=omega;
    S(nanind)=[];
    omegan(nanind)=[];
    s2(n)=spec1d(omegan,S,sqrt(S));
    %[wmax,imax]=max(S);
    %sfold(n)=smooth(s(n),0.05);
    %plot(sfold)
    %[sf,f]=fits(sfold(n),'gaussskew',[100,omegan(imax),0.04,0.01,0],[1 1 1 1 0]);
    %e(n)=f.pvals(2);
    plot(s2(n))
    %pause(0.1)
end
mapplot(s2,h,0.01)
caxis([0,100])
line([0,5],[0 0],'Color','r')