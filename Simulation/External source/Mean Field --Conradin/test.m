q=[0:0.05:1]'*[1 0 0];
%delta=zeros(3,3,4,4,size(q,1));
d=zeros(3,3,4,4,size(q,1));
no=zeros(size(q,1),1,3);
for n=1:size(q,1)
 %   delta(:,:,:,:,n)=dipole_direct(q(n,:),12)-dipole(q(n,:),12);
    d(:,:,:,:,n)=dipole_direct(q(n,:),12);
    no(n,3)=real(sum(d(3,3,:,1,n),3));
    no(n,1)=real(sum(d(1,1,:,1,n),3));
    no(n,2)=real(sum(d(2,2,:,1,n),3));
end

plot(squeeze(q(:,1)),squeeze(no))