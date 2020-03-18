function [pt,ph]=phaseng()
Hcf=cf();
pt=[0:0.05:0.5];
ph=[];
h=[3:0.01:4.8];
for nt=1:length(pt)
    jz=3;
    jx=2;
    jy=0;
    jxyz=[];
    for nh=1:length(h)
        if jz>0.04
            jz=jz-0.01;
        end
        [jx,jy,jz]=remf(h(nh),pt(nt),jx,jy,jz,Hcf);
        jxyz=[jxyz;[jx,jy,jz]];
        h(nh)
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
   
%save pg1 pt ph

load pg1 pt ph

epsilon=0.01;
omega=[0.12:0.02:0.28];
H=[];
omgs=[];
for nh=1:length(ph)
    S=chiqwn([2.1 0 0], omega, epsilon, ph(nh), pt(nh));
    plot(omega,S)
    %waitforbuttonpress()
    H=[H,ph(nh)];
    [m,ind]=max(S);
    omgs=[omgs,omega(ind)];
    nh
end
plot(pt,omgs,'.')