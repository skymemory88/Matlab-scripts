function [fff, freq_total, chi0, chiq] = linear_response(ion,eee,fff,freq_total,ttt,vvv,gama)
% Calculation of susceptibilities
E = eee;
V = vvv;
fields = fff;
% fields = vecnorm(fff);


%Initiate ionJ operators
ionJ = ion.J;
Jz = diag(ionJ:-1:-ionJ); % Jz = -J, -J+1,...,J-1,J
Jp = diag(sqrt((ionJ-((ionJ-1):-1:-ionJ) ).*(ionJ+1+( (ionJ-1):-1:-ionJ) )),1); % electronic spin ladder operator
Jm = Jp'; % electronic spin ladder operator

if ion.hyp > 0 % hyperfine interaction option
    %Initiate I operators
    ionI = ion.I;
    Iz = diag(ionI:-1:-ionI); %Iz = -I, -I+1,...,I-1,I
    IhT.z = kron(eye(2*ionJ+1),Iz); % Expand Hilbert space
    Ip = diag(sqrt((ionI-((ionI-1):-1:-ionI)).*(ionI+1+((ionI-1):-1:-ionI))),1); % Nuclear spin ladder operator
    Im = Ip'; % Nuclear spin ladder operator
    Iph = kron(eye(2*ionJ+1),Ip); % Expand to match the dimension of Hilbert space
    Imh = kron(eye(2*ionJ+1),Im);
    IhT.x = (Iph+Imh)/2;
    IhT.y = (Iph-Imh)/2i;

    JhT.z = kron(Jz,eye(2*ionI+1)); % Expand Jz space to include nuclear degree of freedom
    Jph = kron(Jp,eye(2*ionI+1)); % Expand to match the dimension of Hilbert space
    Jmh = kron(Jm,eye(2*ionI+1));
    JhT.x = (Jph+Jmh)/2;
    JhT.y = (Jph-Jmh)/2i;

    % IJ_hT.x = JhT.x + NUCf/ELEf*IhT.x; % Hybridized electronuclear spin operator
    % IJ_hT.y = JhT.y + NUCf/ELEf*IhT.y;
    % IJ_hT.z = JhT.z + NUCf/ELEf*IhT.z;
else
    JhT.x = (Jp+Jm)/2;
    JhT.y = (Jp-Jm)/2i;
    JhT.z = Jz;

    IhT.x = 0;
    IhT.y = 0;
    IhT.z = 0;
end

chi0 = zeros(3,3,length(freq_total(1,:)),size(fields,2));
for ii = 1:size(freq_total,2) %calculate susceptibility for all frequencies
    freq = freq_total(ii);
    f2E = 1/241.8;  % GHz to meV
    omega = freq*f2E;   % define frequency sweep range (meV)
    parfor jj = 1:size(fields,2) % calculate susceptibility for all fields
        v = squeeze ( squeeze(V(jj,:,:,:)) ); % Obtain the corresponding eigen vectors
        en = squeeze ( squeeze(E(jj,:,:)) ); % Obtain the corresponding eigen energies in meV
%       N = length(en);
%       chi_t = zeros(1,N^2);
%       ll = 1;
        beta = 1/(ttt*8.617E-2); %[meV^-1]
        z = sum(exp(-beta*en));
        zz = exp(-beta*en)/z;
        [n,np] = meshgrid(zz,zz);
        NN = n-np;
        [ee,eep] = meshgrid(en,en);
        EE = eep-ee-omega;
        gamma = ones(size(EE))*gama;
        Gi = gamma ./ (EE.^2 + gamma.^2); 
        Gr = EE ./ (EE.^2 + gamma.^2);  

        ELEf = ion.gLande * 0.05788;     % Lande factor x Bohr magneton (meV T^-1)
        NUCf = 4.732 * 3.1519e-5;   % Nuclear Lande factor
        JxhT = JhT.x * ELEf;
        IxhT = IhT.x * NUCf;

        JyhT = JhT.y * ELEf;
        IyhT = IhT.y * NUCf;

        JzhT = JhT.z * ELEf;
        IzhT = IhT.z * NUCf;

        ttx  = v'  * (JxhT+IxhT) * v; 
        tty  = v'  * (JyhT+IyhT) * v; 
        ttz  = v'  * (JzhT+IzhT) * v;

% Calculate susceptibilities xx
        chii  = (ttx) .* (ttx.') .* NN .* Gi;
        chir = (ttx) .* (ttx.') .* NN .* Gr;
        xxi =  real( sum(sum(chii)) )   ;
        xxr =  real( sum(sum(chir)) )  ; 

% Calculate susceptibilities xy
        chii  = (ttx) .* (tty.') .* NN .* Gi;
        chir = (ttx) .* (tty.') .* NN .* Gr;
        xyi =  real( sum(sum(chii)) )   ;
        xyr =  real( sum(sum(chir)) )  ;   

% Calculate susceptibilities xz
        chii  = (ttx) .* (ttz.') .* NN .* Gi;
        chir = (ttx) .* (ttz.') .* NN .* Gr;
        xzi =  real( sum(sum(chii)) )   ;
        xzr =  real( sum(sum(chir)) )  ;

% Calculate susceptibilities yx
        chii  = (tty) .* (ttx.') .* NN .* Gi;
        chir = (tty) .* (ttx.') .* NN .* Gr;
        yxi =  real( sum(sum(chii)) )   ;
        yxr =  real( sum(sum(chir)) )  ; 

% Calculate susceptibilities yy
        chii  = (tty) .* (tty.') .* NN .* Gi;
        chir = (tty) .* (tty.') .* NN .* Gr;
        yyi =  real( sum(sum(chii)) )   ;
        yyr =  real( sum(sum(chir)) )  ;   

% Calculate susceptibilities yz
        chii  = (tty) .* (ttz.') .* NN .* Gi;
        chir = (tty) .* (ttz.') .* NN .* Gr;
        yzi =  real( sum(sum(chii)) )   ;
        yzr =  real( sum(sum(chir)) )  ;

% Calculate susceptibilities zx
        chii  = (ttz) .* (ttx.') .* NN .* Gi;
        chir = (ttz) .* (ttx.') .* NN .* Gr;
        zxi =  real( sum(sum(chii)) )   ;
        zxr =  real( sum(sum(chir)) )  ; 

% Calculate susceptibilities zy
        chii  = (ttz) .* (tty.') .* NN .* Gi;
        chir = (ttz) .* (tty.') .* NN .* Gr;
        zyi =  real( sum(sum(chii)) )   ;
        zyr =  real( sum(sum(chir)) )  ;   

% Calculate susceptibilities zz
        chii  = (ttz) .* (ttz.') .* NN .* Gi;
        chir = (ttz) .* (ttz.') .* NN .* Gr;
        zzi =  real( sum(sum(chii)) )   ;
        zzr =  real( sum(sum(chir)) )  ;

    xr = [xxr xyr xzr
          yxr yyr yzr
          zxr zyr zzr]; % Real part of susceptibility
    xi = [xxi xyi xzi
          yxi yyi yzi
          zxi zyi zzi]; % Imaginary part of susceptibility

    chi0(:,:,ii,jj) = xr + 1i.*xi;
    end
end
[~, ~, chiq] = RPA([0 0 0], fff, freq_total, ion, chi0); % calculate RPA susceptibility
end

function [continu_var, freq_total, chiq] = RPA(qvec, continu_var, freq_total, ion, chi0)
unitN = 1; % Number of magnetic atoms in unit cell
% muB = 9.274e-24; % [J/T]
% mu0 = 4e-7*pi; % [H/m]
% J2meV = 6.24151e+21; % [mev/J]
% gfac = ion.gLande^2*muB^2*mu0/4/pi*J2meV; % (gL * const.muB)^2 * const.mu0/(4pi) [meV.m^3]

chiq = zeros(3,3,length(freq_total(1,:)),size(continu_var,2),size(qvec,1));
D = zeros(3,3,unitN,unitN,size(qvec,1));
for jj = 1:size(qvec,1)
%     D(:,:,:,:,jj) = gfac*10^30*dipole_direct(qvec(jj,:),const.dpRng,lattice)...
%         + exchange(qvec(jj,:),ion.ex(2),lattice,tau); % [meV]
    D(:,:,:,:,jj) = exchange(qvec(jj,:), ion.ex, ion.abc, ion.tau); % [meV]
end

deno = zeros(3,3,size(freq_total,1),size(continu_var,2)); % RPA correction factor (denominator)
for ii = 1:size(continu_var,2)
    for nq = 1:size(qvec,1)
        Jq = sum(sum(D(:,:,:,:,nq),4),3)/unitN; % average over the all sites in the unit cell and avoid double counting [meV]
        parfor kk = 1:length(freq_total(1,:))
            MM = chi0(:,:,kk,ii).*Jq; %[meV*meV^-1]
            deno(:,:,kk,ii) = (ones(size(MM))- MM); % Suppress parfor warning for this line
            chiq(:,:,kk,ii,nq) = squeeze(deno(:,:,kk,ii)).\chi0(:,:,kk,ii);
        end
    end
end

% subs = ["xx", "yy", "zz"];
% for ii = 1:3
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze((Mx(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     
%     figure
%     hp0 = pcolor(continu_var,freq_total,squeeze(abs(deno(ii,ii,:,:))));
%     set(hp0, 'edgeColor','none')
%     xlabel('Magnetic field (T)')
%     ylabel('Frequency (GHz)')
%     for nq = 1:size(qvec,1)
%         figure
%         hp0 = pcolor(continu_var,freq_total,real(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Real part of \chi_',char(subs(ii)),'^{RPA}'])
%         figure
%         hp0 = pcolor(continu_var,freq_total,imag(squeeze(chiq(ii,ii,:,:,nq))));
%         set(hp0, 'edgeColor','none');
%         colorbar
%         xlabel('Magnetic field (T)')
%         ylabel('Frequency (GHz)')
%         title(['Imaginary part of \chi_',char(subs(ii)),'^{RPA}'])
%     end
% end
end

function d = exchange(q, Jex, a, tau, nn)

% This function performs a brute force summation of
% the q-dependent exchange coupling fo a non-Bravais lattice.
% q=[h k l] is the q-vector given in Miller-indicies.
% nr is the number of unit cells that should be summed
% in each direction.

% For N moments in the unit cell, there will be N 
% coupling parameters J_ij. Many of these will be symmetry related
% (e.g. J_ij=J_ji), so we just calculate J_1j. 
% The result is a (3x3xN) matrix, where the first two dimensions
% are the x,y and z components. The last dimension holds the coupling
% between different ions in the unit cell.
switch nargin
    case 4    
        N = 1;
    case 5
        N = nn;
    otherwise
        error('Incorrect number of input argument for exchange()!')
end

% Convert tau to a
tau = tau*a;

% Unit cell volume
vol = sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b = [2*pi*cross(a(2,:),a(3,:))/vol
     2*pi*cross(a(3,:),a(1,:))/vol
     2*pi*cross(a(1,:),a(2,:))/vol];
% Convert q from Miller indicies to reciprocal angstroms
q = q*b;
% % Length of q
% qq=sqrt(sum(q.*q));

% % Kronecker delta in x,y,and z
% delta=[1 0 0;0 1 0;0 0 1]; % replaced by built-in eq(m,n) function, -Yikai 2021.03.08

[x,y,z] = meshgrid(-N:N,-N:N,-N:N);
hkl = [x(:) y(:) z(:)];

d = zeros(3,3,size(tau,1),size(tau,1));
for nt = 1:size(tau,1)
    for mt = 1:nt
        r = hkl*a;
        r(:,1) = r(:,1)-tau(nt,1)+tau(mt,1);
        r(:,2) = r(:,2)-tau(nt,2)+tau(mt,2);
        r(:,3) = r(:,3)-tau(nt,3)+tau(mt,3);
        rr = sum(r.^2,2);
        r(rr > (N * 5.175)^2 | rr<0.01,:) = []; % Remove points outside of cut-off range and singularities
        rr(rr > (N * 5.175)^2 | rr<0.01) = [];
        exp_qr = exp(-1i*q*r');
        for n = 1:3
            for m = 1:3
                d(n,m,nt,mt) = exp_qr*(rr<(N*10.75)^2)*Jex*eq(n,m); % exchange interaction
%                 d(n,m,nt,mt) = exp_qr*(rr<14)*Jex*delta(n,m); % original code, unsure of the rr<14
            end
        end
        %  d(:,:,nt,mt) = d(:,:,nt,mt)+(4*pi/3)*0.01389*eye(3)/4; %Lorentz
        d(:,:,mt,nt) = conj(d(:,:,nt,mt));
    end
end
% Sometimes the coupling is given in units of the unit cell volume:
%d=d*vol;

return
end