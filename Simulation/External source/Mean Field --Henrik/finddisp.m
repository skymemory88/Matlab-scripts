%load chiq_results_save

wdisp=zeros(length(qz),length(h));
for n=1:length(qz)
  for m=1:length(h)
    [a,b]=sort(diff(squeeze(real(chiqz(1,1,:,n,m)))));
    wdisp(n,m)=(w(b(2))+w(b(2)+1))/2;
  end
end    

q=qz;
chi=chiqz;
chim=chiqzm;
Skw=zeros(length(w),length(qz),length(h));
for nq=1:length(q)
  for nw=1:length(w)
    for nh=1:length(h)
  Skw(nw,nq,nh)=1/pi* ...
      sum(sum((eye(3)-(q(nq,:)'*q(nq,:))/(q(nq,:)*q(nq,:)')).* ...
              imag(chi(:,:,nw,nq,nh)-chim(:,:,nw,nq,nh))));
end
end
end

wdisp=zeros(length(qz),length(h));
for n=1:length(qz)
  for m=1:length(h)
    [a,b]=sort(diff(squeeze(real(chiqz(1,1,:,n,m)))));
    wdisp(n,m)=(w(b(2))+w(b(2)+1))/2;
  end
end    
