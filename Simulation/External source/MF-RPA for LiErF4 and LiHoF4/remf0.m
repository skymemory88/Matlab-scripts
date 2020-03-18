function [momente,e]=remf(h,t,momente,Jex,hyper)


if nargin<4
    Jex=0;
end

if size(h,2)==1
    hvec=[0 0 h];
else
    hvec=h;
end


J=15/2;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

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
else
    hJx=Jx;
    hJy=Jy;
    hJz=Jz;
end


%Lorenzterm
Lorenz=4.8813/1000;
eins=zeros(3,3,4,4); eins(1,1,:,:)=1; eins(2,2,:,:)=1; eins(3,3,:,:)=1;

%d=(6/5)^2*0.05368*(dipole_direct([0 0 0],15)+exchange([0,0,0],Jex));
d=(6/5)^2*0.05368*(dipole_direct([0 0 0],15)+exchange([0,0,0],Jex))+eins*Lorenz/4;

for iterations=1:200
    % iterations
    % Diagonalize
    momente_old=momente;
%     for ionn=4:-1:1
    for ionn=1:4
        h_dipol=[0 0 0];
        for ionm=1:4
            h_dipol=h_dipol+momente_old(ionm,:)*d(:,:,ionm,ionn)';
        end
        Ham=REHam(hvec,h_dipol,hyper);
        [v,e]=eig(Ham);
        e=real(diag(e));
        e=e-min(e);
        [e,n]=sort(e);
        v=v(:,n);
        
        % e contains eigenvalues,
        % rows of v contains eigenvectors
        if t==0
            % At zero temperature, use only lowest eigenvalue.
            jx=real(v(:,1)'*hJx*v(:,1));
            jy=real(v(:,1)'*hJy*v(:,1));
            jz=real(v(:,1)'*hJz*v(:,1));
        else
            % Boltzman factor (with t in Kelvin)
            % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
            e=e-min(e);
            z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));
            jx=real(diag(v'*hJx*v))'*z;
            jy=real(diag(v'*hJy*v))'*z;
            jz=real(diag(v'*hJz*v))'*z;
        end
        momente(ionn,:)=[jx,jy,jz];
    end
    % Itterate untill absolute change is less than 0.001
    
    if sum(sum(abs(momente_old-momente)))<0.0001
%        iterations
        break
    end
    
end
iterations
return



