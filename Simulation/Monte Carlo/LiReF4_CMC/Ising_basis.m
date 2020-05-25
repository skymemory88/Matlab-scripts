function [VV,H]=Ising_basis(field, ion, Jx, Jy, Jz)
% computes the effective Ising basis for the given external magnetic field

muB=0.05788; % Bohr magneton, in meV
gLande=ion.gLande; % Landé factor

if(nargin<3)
    J=ion.J;
    Jz=diag(J:-1:-J);
    Jp=diag(sqrt((J-[J-1:-1:-J]).*(J+1+[J-1:-1:-J])),1);
    Jm=Jp';
    Jx=(Jp+Jm)/2;
    Jy=(Jp-Jm)/2i;
end

%diagonalization of Hcf+Hzee
H=ion.Hcf-muB*gLande*(field(1)*Jx+field(2)*Jy+field(3)*Jz);
[V,~]=eig(0.5*(H+H'));

% %diagonalization of Hcf+Hzee and sort by eigen-energy
% H=ion.Hcf-muB*gLande*(field(1)*Jx+field(2)*Jy+field(3)*Jz);
% [V,e]=eig(0.5*(H+H'));
% [~,idx] = sort(e);
% clearvars e
% V = V(:,idx);

% diagonalization of the projection of Jz on the subspace
[Vp,~]=eig(V(:,1:2)'*Jz*V(:,1:2));

% Ising basis
VV=V(:,1:2)*[Vp(:,1) Vp(:,2)];

end