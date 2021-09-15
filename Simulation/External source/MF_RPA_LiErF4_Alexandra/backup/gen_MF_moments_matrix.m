function [mJx,mJy,mJz,e,z,v] = gen_MF_moments_matrix(ion,hvec,h_mf,t,ishf)

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

%% Calculates matrix expected value
  mJx=v'*JIxh*v;
  mJy=v'*JIyh*v;
  mJz=v'*JIzh*v;

if t==0
    z=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0  
    % At zero temperature, use only lowest eigenvalue.
else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e*beta)/sum(exp(-e*beta)); % Density matrix element
end

