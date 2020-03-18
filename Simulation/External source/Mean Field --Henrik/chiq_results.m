clear all
w=0:0.01:0.5;
h=0:1:6;
%qx=[0 0.01 0.05:0.05:0.95 0.99 1]'*[1 0 0];
%qy=[0 0.01 0.05:0.05:0.95 0.99 1]'*[0 1 0];
qz=[0 0.01 0.05:0.05:0.95 0.99 1 1.01 1.05:0.05:1.95 1.99 2]'*[0 0 1];

%chiqx=zeros(3,3,length(w),size(qx,1),length(h));
%chiqy=zeros(3,3,length(w),size(qy,1),length(h));
chiqz=zeros(3,3,length(w),size(qz,1),length(h));

t=0;
epsilon=0.01;

for n=1:length(h)
  h(n)
%  chiqx(:,:,:,:,n)=chi_qw(qx,w,epsilon,h(n),t);
%  chiqy(:,:,:,:,n)=chi_qw(qy,w,epsilon,h(n),t);
  chiqz(:,:,:,:,n)=chi_qw(qz,w,epsilon,h(n),t);
%  chiqxm(:,:,:,:,n)=chi_qw(qx,-w,epsilon,h(n),t);
%  chiqym(:,:,:,:,n)=chi_qw(qy,-w,epsilon,h(n),t);
  chiqzm(:,:,:,:,n)=chi_qw(qz,-w,epsilon,h(n),t);
  save chiq_results_save2
end

