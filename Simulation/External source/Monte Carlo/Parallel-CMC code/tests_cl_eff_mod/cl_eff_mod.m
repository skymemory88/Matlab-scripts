function [mom, E, params]=cl_eff_mod(alpha,beta,field,ion,params)
% computes the effective moment and energy for an ion, using the
% classical model derived for the quantum effective model ie computes the
% mean value of the effective quantum momentum and hamiltonian in the state
% defined by alpha and beta

% muB=0.05788; % Bohr magneton, in meV/T
J=ion.J; % Spin
% gLande=ion.gLande; % Landé factor

% definition of the J operators
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

% % crystal field hamiltonian
% Hcf=ion.Hcf;
% 
% % diagonalization of Hcf+Hzee
% H=Hcf-muB*gLande*(field(1)*Jx+field(2)*Jy+field(3)*Jz);
% 
% [V,~]=eig(0.5*(H+H'));
% 
% % diagonalization of the projection of Jz on the subspace
% [Vp,~]=eig(V(:,[1 2])'*Jz*V(:,[1 2]));
%  
% % projection and computation of J and H in the new basis
% VV=V(:,1:2)*[Vp(:,1) Vp(:,2)];


if(params.field_changed==1)
    [VV,H]=Ising_basis(field, ion, Jx, Jy, Jz);
    params.field_changed=0;
else
    VV=ion.VV;
    muB=0.05788; % Bohr magneton, in meV
    H=ion.Hcf-muB*ion.gLande*(field(1)*Jx+field(2)*Jy+field(3)*Jz);
end

Jxt=VV'*Jx*VV;
Jyt=VV'*Jy*VV;
Jzt=VV'*Jz*VV;
Hzeff=VV'*H*VV;

% computation of the mean values on the state defined by alpha and beta, we
% take alpha/2 to avoid an asymmetry between +z and -z.
state=[cos(alpha/2) ; sin(alpha/2).*exp(1i*beta)];

mom(1)=real(state'*Jxt*state);
mom(2)=real(state'*Jyt*state);
mom(3)=real(state'*Jzt*state);
E=real(state'*Hzeff*state);

end
