
for n=1:length(h)
  [Ham,jx,jy,jz,e(:,n),v]=meanfield(h(n),0,0,0,1);
  [split,vg,ve,Hcfz,v,ecf(:,n)]=crystfield(h(n));
end

h1=plot(h,ecf(1:3,:))
hold on
  h2=plot(h,e(1:3,:),'--')
hold off
nicefig
axis([0 6 -1 2])
xlabel('Transverse field [T]')
ylabel('Energy [meV]')
title('Crystal field levels in LiHoF_4')
legend([h1(1) h2(1)],'V_{cf} only','V_{cf} and J_{mf}',2)
print -depsc crystfield_levels.eps

h1=plot(h,ecf(2,:)-ecf(1,:),h,ecf(3,:)-ecf(1,:))
hold on
  h2=plot(h,e(2,:)-e(1,:),'--',h,e(3,:)-e(1,:),'--')
hold off
nicefig
axis([0 6 -0.1 1.85])
xlabel('Transverse field [T]')
ylabel('Energy [meV]')
title('Crystal field splitting in LiHoF_4')
legend([h1(1) h2(1)],'V_{cf} only','V_{cf} and J_{mf}',2)
print -depsc crystfield_split.eps




