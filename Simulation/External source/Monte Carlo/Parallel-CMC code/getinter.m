function inter=getinter(params)

%Lorenzterm = N/V * mu_0 * (mu_b)^2/4pi * 4pi/3
% Lorentz=3.124032/(1000); %Ho
% 0.000745836=(N/V)*0.05368

t=tic;

L=params.L;

[n,m,l,k,h]=ndgrid(1:4,1:4,-L:L,-L:L,-L:L); % Generate all possible combinations of ion pairs
hklmn_d=[h(:) k(:) l(:) m(:) n(:)];

inter=zeros(3,3,size(hklmn_d,1));

% % Er in Unit cell volume
% V_uc=params.abc(1,1)*params.abc(2,2)*params.abc(3,3); 
% N_tot=(Nx+1)*(Ny+1)*(Nz+1)*4;
% N=4/V_uc; % Number of ions per unit volume
% V=V_uc/4; % Volume per ion
% Lorentz=3.15452/1000; %Er
% demagn=ellipsoid_demagn(1);

% % Lorenzterm = N/V * mu_0 * (mu_b)^2/4pi * 4pi/3
V=params.abc(1,1)*L*params.abc(2,2)*L*params.abc(3,3)*L; % Box volume = lattice constants * # of unit cells
% N_tot=(L^3)*4; % Number of ions
% demagn=(N_tot/V)*ellipsoid_demagn(1);
Lorentz=(1/V)*(4*pi/3); % why is it not (N_tot/V) ? --Yikai
ct=0.05368; % ct = mu_0*mu_B^2/4pi ([meV*Ang^-3])
% Lorentz=1/(params.abc(1,1)*params.replicas(1)*L)^3;

% iterate through the four ions inside each unit cell
for i=1:size(hklmn_d,1)
    r_ij=[h(i);k(i);l(i)]+params.pos(:,n(i))-params.pos(:,m(i)); % vector connecting tow ions
    idx=latmod(params,r_ij);
    % inter(:,:,idx)=dipole_direct(params.abc*r_ij,params.abc*[L; L; L],params.replicas);
     inter(:,:,idx)=ct*(dipole_direct(params.abc*r_ij,params.abc*[L; L; L],params.replicas)+(eye(3)*Lorentz));
%     inter(:,:,idx)=ct*(dipole_direct(params.abc*r_ij,params.abc*[L; L; L],params.replicas)+(eye(3)*Lorentz)-pi*demagn); %Er
%     inter(:,:,idx)=1.2*ct*(dipole_direct(params.abc*r_ij,params.abc*[L; L; L],params.replicas)+(eye(3)*Lorentz)-pi*demagn); %OK for 1 unit cell
end

toc(t);

end