function [mom, E, field_change]=cl_eff_mod(alpha,beta,field,ion,field_change)
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

if (field_change==true)
    [VV,H]=Ising_basis(field, ion, Jx, Jy, Jz);
    field_change=false;
else
    VV=ion.VV;
    muB=0.05788; % Bohr magneton, in meV
    H=ion.Hcf-muB*ion.gLande*(field(1)*Jx+field(2)*Jy+field(3)*Jz); % Crystal field and Zeeman terms
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
