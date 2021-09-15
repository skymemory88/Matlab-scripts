function [jx,jy,jz,e,z,v] = gen_MF_moments(ion,hvec,h_mf,t,ishf)

[Ham,Jxh,Jyh,Jzh,Ixh,Iyh,Izh] =calc_ham(ion,hvec,h_mf,ishf);

%Diagonalize
[v,e]=eig(Ham);
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);
beta = 11.6/t; % [kB*T]^-1 in [meV]
%Calculate Matrix elements and Moments
    JIxh = Jxh;
    JIyh = Jyh;
    JIzh = Jzh;
% end

if t==0
    z=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0  
    % At zero temperature, use only lowest eigenvalue.
    jx=real(v(:,1)'*JIxh*v(:,1));
    jy=real(v(:,1)'*JIyh*v(:,1));
    jz=real(v(:,1)'*JIzh*v(:,1));
else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e*beta)/sum(exp(-e*beta)); % Density matrix element
    jx=real(diag(v'*JIxh*v))'*z;
    jy=real(diag(v'*JIyh*v))'*z;
    jz=real(diag(v'*JIzh*v))'*z;
end
