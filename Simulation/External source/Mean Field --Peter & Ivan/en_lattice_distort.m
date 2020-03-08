fm=[0 0 1
    0 0 1
    0 0 1
    0 0 1]*3.5;

afm=[-1 0 0
    1 0 0
    1 0 0
    -1 0 0]*3.5;

ca=1:0.01:2.5;
J=7/2;
gL=2;
B=[0,0,0.3];
enfm=zeros(size(ca));
enafm=zeros(size(ca));

for r=1:length(ca)
    abc=diag([1, 1, ca(r)]*5.162);
    dipole=dipole_direct([0,0,0],15,abc);
%     Jz=diag(J:-1:-J);
%     Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
%     Jm=Jp';
%     Jx=(Jp+Jm)/2;
%     Jy=(Jp-Jm)/2i;
%     Hdfm=zeros([4,2*J+1,2*J+1]);
%     Hdafm=zeros([4,2*J+1,2*J+1]);
    denfm=[0,0,0,0];
    denafm=[0,0,0,0];
    for n=1:4
        for m=1:4
            dfm=dipole(:,:,n,m)*fm(n,:)';
            denfm(n)=denfm(n)+0.05368*gL*(dfm(1)*fm(m,1)+dfm(2)*fm(m,2)+dfm(3)*fm(m,3));
            dafm=dipole(:,:,n,m)*afm(n,:)';
            denafm(n)=denafm(n)+0.05368*gL*(dafm(1)*afm(m,1)+dafm(2)*afm(m,2)+dafm(3)*afm(m,3));
        end
        denfm(n)=denfm(n)+0.05788*gL*(fm(n,1)*B(1)+fm(n,2)*B(2)+fm(n,3)*B(3));
        denafm(n)=denafm(n)+0.05788*gL*(afm(n,1)*B(1)+afm(n,2)*B(2)+afm(n,3)*B(3));
    end
    enfm(r)=-sum(denfm);
    enafm(r)=-sum(denafm);
end
figure
plot(ca,enfm,ca,enafm)
title('LiGdF_4')
xlabel('c/a');
ylabel('E');
legend('FM','AFM');