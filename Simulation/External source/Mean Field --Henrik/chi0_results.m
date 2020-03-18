w=0:0.05:2;
h=0:0.5:6;

chi0mf=zeros(3,3,length(w),length(h));
chi0cf=zeros(3,3,length(w),length(h));

t=0;
epsilon=0.01;

for n=1:length(w)
for m=1:length(h)
  chi0mf(:,:,n,m)=chi0_w(h(m),t,w(n),epsilon,1);
  chi0cf(:,:,n,m)=chi0_w(h(m),t,w(n),epsilon,0);
end
end