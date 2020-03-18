function Ham=REHam(hvec,h_dipol,hyper)
%H0 bestimmen

% %Initiate J operators
J=15/2;
Jid=eye(2*J+1);
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

Hcf=cf_Er;

if hyper~=0
    I=7/2;
    Iid=eye(2*I+1);
    Iz=diag(I:-1:-I);
    Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
    Im=Ip';
    Ix=(Ip+Im)/2;
    Iy=(Ip-Im)/2i;
    hJx=kron(Jx,Iid);
    hJy=kron(Jy,Iid);
    hJz=kron(Jz,Iid);
    hIx=kron(Jid,Ix);
    hIy=kron(Jid,Iy);
    hIz=kron(Jid,Iz);
    hHcf=kron(Hcf,Iid);
    
    Hzeeman=(-1.2*0.05788)*(hvec(1)*hJx+hvec(2)*hJy+hvec(3)*hJz);
    
    Ham=hHcf+Hzeeman-h_dipol(1)*hJx-h_dipol(2)*hJy-h_dipol(3)*hJz +hyper*(hJx*hIx+hJy*hIy+hJz*hIz);
    
else
    
    Hzeeman=(-1.2*0.05788)*(hvec(1)*Jx+hvec(2)*Jy+hvec(3)*Jz);
    
    Ham=Hcf+Hzeeman-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;
end



