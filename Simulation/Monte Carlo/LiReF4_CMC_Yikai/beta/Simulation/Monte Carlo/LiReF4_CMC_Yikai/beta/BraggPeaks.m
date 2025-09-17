function BraggPeaks(pos, spin, wavelength, kVec, Xaxis, xrange)
% Simulate neutron diffraction from a REAL finite lattice given by pos/spin.
% No periodic replication. Bragg-like peaks emerge from the finite-lattice
% structure factor (shape-transformed).
%
% pos:  N x 3 (Angstrom)
% spin: N x 3 (magnetic moments, e.g. in mu_B, can be zero rows for nuclear-only)
% wavelength: neutron lambda [Angstrom]
% kVec: incident wavevector (1/Angstrom). [] => POWDER; nonempty => SINGLE-CRYSTAL
% Xaxis: 'angle' (2theta) or 'space' (d-spacing). (Only affects plotting axis.)
% xrange: grid of 2theta (deg) for 'angle' or d (Angstrom) for 'space'.
%
% Notes:
% - Uses element-independent b_j (set scalar below); supply a vector b_j if you have it.
% - Magnetic form factor is set to Er3+ as an example; replace as needed.

% ---- basic checks ----
if nargin < 6, xrange = []; end
if isempty(spin), spin = zeros(size(pos)); end
pos  = reshape(pos,  [], 3);
spin = reshape(spin, [], 3);
if size(pos,1) ~= size(spin,1), error('pos and spin must have same N.'); end
if ~isempty(kVec) && numel(kVec) ~= 3, error('kVec must be 1x3 or [].'); end
if ~ismember(Xaxis, {'space','angle'}), error('Xaxis must be ''space'' or ''angle''.'); end

% ---- physical constants & options ----
r0       = 2.818e-15;          % m, classical electron radius
gamma_n  = -1.913;             % neutron gyromagnetic ratio (dimensionless)
Cmag     = abs(gamma_n)*r0/2;  % m, magnetic scattering length prefactor
b_nucl   = 7.0e-15;            % m, coherent nuclear b (set per atom vector for realism)
n_pwdir  = 256;                % # directions for magnetic powder averaging
n_phi    = 256;                % # azimuth samples on the Ewald ring (single crystal)
% You can increase n_pwdir/n_phi for sharper peaks (cost ~linear).

% ---- build Q-grid according to axis choice ----
if strcmpi(Xaxis,'angle')
    if isempty(xrange), twoTheta = linspace(5,120,1200); else, twoTheta = xrange(:)'; end
    Q = (4*pi/wavelength) * sind(twoTheta/2);  % 1/Angstrom
    x_for_plot = twoTheta;
    xlab = '2\theta (deg)';
else % 'space' -> d-spacing axis
    if isempty(xrange), d = linspace(0.7,6,2000); else, d = xrange(:)'; end
    d(d==0) = eps;
    Q = 2*pi ./ d;                              % 1/Angstrom
    x_for_plot = d;
    xlab = 'd (\AA)';
end

% ---- choose mode: powder if kVec empty, else single-crystal ----
isPowder = isempty(kVec);

if isPowder
    % ======================== POWDER (real finite lattice) =======================
    % Nuclear: exact Debye powder average (NO replication, preserves peak positions)
    % Magnetic: MC average over Q-direction (robust with arbitrary spin textures)
    [I_N, I_M] = powder_intensity_real(pos, spin, Q, b_nucl, Cmag, n_pwdir, @magnetic_form_factor);
    I_total = I_N + I_M;

    figure; semilogy(x_for_plot, max(I_total,1e-30),'g-','LineWidth',1.4); hold on;
    semilogy(x_for_plot, max(I_N,1e-30), 'b--','LineWidth',1.1);
    semilogy(x_for_plot, max(I_M,1e-30), 'r--','LineWidth',1.1);
    xlabel(xlab,'Interpreter','tex'); ylabel('Intensity (arb, barns scale)');
    title('Powder diffraction from finite real lattice (nuclear + magnetic)');
    legend('Total','Nuclear','Magnetic','Location','best'); grid on; hold off;

else
    % ===================== SINGLE-CRYSTAL (Ewald ring sampling) ==================
    % --- single-crystal setup (make |k| consistent with wavelength) ---
    khat = kVec(:).' / norm(kVec);          % direction from input
    k_i  = 2*pi / wavelength;               % magnitude from lambda

    % orthonormal basis perpendicular to k_i
    U = null(khat);  % 3x2
    e1 = U(:,1)'; e2 = U(:,2)';

    [I_N, I_M] = single_intensity_real(pos, spin, Q, b_nucl, Cmag, k_i, khat, e1, e2, n_phi, @magnetic_form_factor);
    I_total = I_N + I_M;

    figure; semilogy(x_for_plot, max(I_total,1e-30),'g-','LineWidth',1.4); hold on;
    semilogy(x_for_plot, max(I_N,1e-30), 'b--','LineWidth',1.1);
    semilogy(x_for_plot, max(I_M,1e-30), 'r--','LineWidth',1.1);
    xlabel(xlab,'Interpreter','tex'); ylabel('Intensity (arb, barns scale)');
    title('Single-crystal diffraction (Ewald sphere) from finite real lattice');
    legend('Total','Nuclear','Magnetic','Location','best'); grid on; hold off;
end
end

% ============================== POWDER CORE =====================================
function [I_N, I_M] = powder_intensity_real(pos, spin, Q, b_nucl, Cmag, n_pwdir, fmag_handle)
% Nuclear: exact Debye formula (pair sum) -> sharp powder peaks from your finite lattice.
% Magnetic: sample Q-hat directions and average |F_M|^2 (exact for given N, M).

NQ = numel(Q);
N  = size(pos,1);
I_N = zeros(NQ,1);
I_M = zeros(NQ,1);

% --- Debye nuclear (exact powder avg) ---
% I_N(Q) = sum_j b_j^2 + 2*sum_{j<k} b_j b_k * sin(Q r_jk)/(Q r_jk)
% Here b_j is scalar b_nucl for all; change to per-atom vector if you have it.
rij = pdist(pos);                          % Angstrom
bj2 = N * (b_nucl^2);                      % sum_j b_j^2 (all equal)
% vectorized pair contribution weights:
pair_w = (b_nucl^2) * 2;                   % 2*b_j*b_k for all equal

for i = 1:NQ
    Qi = Q(i);
    if Qi == 0
        I_N(i) = bj2;  %#ok<*AGROW>
    else
        x = Qi * rij;
        I_N(i) = bj2 + pair_w * sum(sin(x)./x);
    end
end

% Convert m^2 -> barns (for readability)
I_N = I_N * 1e28;

% --- Magnetic (MC over Q-hat) ---
if any(spin(:) ~= 0)
    for i = 1:NQ
        Qi = Q(i);
        if Qi == 0, I_M(i) = 0; continue; end
        % average over random directions of Q-hat:
        acc = 0;
        for m = 1:n_pwdir
            % random unit vector qhat (uniform on sphere)
            v = randn(1,3); v = v / norm(v);
            qvec = Qi * v;
            phase = exp(1i * (pos * qvec.'));        % Nx1
            % project spins perpendicular to Q
            spin_perp = spin - (spin * v.') .* v;    % Nx3
            % magnetic form factor (can be ion-specific)
            fQ = fmag_handle(Qi);
            Fm_vec = Cmag * fQ * sum(phase .* spin_perp, 1);   % 1x3
            acc = acc + sum(abs(Fm_vec).^2);
        end
        I_M(i) = acc / n_pwdir;
    end
    I_M = I_M * 1e28;  % m^2 -> barns
else
    I_M(:) = 0;
end
end

% =========================== SINGLE-CRYSTAL CORE ================================
function [I_N, I_M] = single_intensity_real(pos, spin, Q, b_nucl, Cmag, k_i, khat, e1, e2, n_phi, fmag_handle)
% For each |Q|, sample the elastic Ewald ring azimuth phi exactly:
% cos(psi) = -Q/(2k), qhat = cospsi*khat + sinpsi*(cosphi*e1 + sinphi*e2)

NQ = numel(Q);
I_N = zeros(NQ,1);
I_M = zeros(NQ,1);

for i = 1:NQ
    Qi = Q(i);
    if Qi <= 0 || Qi > 2*k_i
        I_N(i) = 0; I_M(i) = 0; continue;
    end
    cospsi = -Qi/(2*k_i);
    if abs(cospsi) > 1, continue; end
    sinpsi = sqrt(max(0,1 - cospsi^2));

    Fn_acc = 0; Fm_acc = 0;
    for j = 1:n_phi
        phi  = 2*pi*(j-1)/n_phi;
        qhat = cospsi*khat + sinpsi*(cos(phi)*e1 + sin(phi)*e2);
        qvec = Qi * qhat;

        phase = exp(1i * (pos * qvec.'));          % Nx1

        % Nuclear amplitude and intensity
        Fn  = b_nucl * sum(phase);
        Fn_acc = Fn_acc + abs(Fn)^2;

        % Magnetic
        if any(spin(:) ~= 0)
            spin_perp = spin - (spin * qhat.') .* qhat;  % Nx3
            fQ = fmag_handle(Qi);
            Fm_vec = Cmag * fQ * sum(phase .* spin_perp, 1);  % 1x3
            Fm_acc = Fm_acc + sum(abs(Fm_vec).^2);
        end
    end

    I_N(i) = Fn_acc / n_phi;
    I_M(i) = Fm_acc / n_phi;
end

% m^2 -> barns
I_N = I_N * 1e28;
I_M = I_M * 1e28;
end

% ======================== Magnetic form factor example ==========================
function f_mag = magnetic_form_factor(q)
% Example: Er3+ (4f^11) analytical fit: f(q) = A*exp(-a s) + B*exp(-b s) + C, with s = (q/4π)^2
s = (q/(4*pi))^2;
A = 0.6894; a = 12.12;
B = 0.4115; b = 4.838;
C = -0.1359;
f_mag = A*exp(-a*s) + B*exp(-b*s) + C;
f_mag = max(f_mag, 1e-3);   % clip to avoid negatives at large Q
end
