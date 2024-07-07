function nSpin = therm_nuc_GPU(const, ion, params, beta, field, eSpin, nSpin)
numIso = size(params.isoIdx,1);

% Extract components for all isotopes
spx = eSpin(params.isoIdx, 1);
spy = eSpin(params.isoIdx, 2);
spz = eSpin(params.isoIdx, 3);

spxi = nSpin(params.isoIdx, 1);
spyi = nSpin(params.isoIdx, 2);
spzi = nSpin(params.isoIdx, 3);

Ix = gpuArray(ion.Ix);
Iy = gpuArray(ion.Iy);
Iz = gpuArray(ion.Iz);

% Compute energy contributions for all isotopes at once
EzI = -ion.nLande(ion.idx) * const.muN * const.J2meV * (field(1) * spxi + field(2) * spyi + field(3) * spzi);
E_hyp = ion.A(ion.idx) * (spx .* spxi + spy .* spyi + spz .* spzi);
E_old = EzI + E_hyp;

% Define Hamiltonians (assuming constants are scalars or compatible matrices)
H_hyp =  arrayfun(@(ii) (spx(ii) * Ix + spy(ii) * Iy + spz(ii) * Iz), 1:numIso, 'UniformOutput', false);
H_hyp = ion.A(ion.idx) * cat(3,H_hyp{:});
HzI = -ion.nLande(ion.idx) * const.muN * const.J2meV * (field(1) * Ix + field(2) * Iy + field(3) * Iz);
ham_nuc = H_hyp + HzI;

% Initialize GPU arrays for new spins and energies
newSpins = gpuArray.zeros(numIso, 3);
E_new = gpuArray.zeros(numIso, 1);
prob = gpuArray.zeros(numIso, 1);

% Loop over isotopes for eigen decomposition and random configuration (not vectorized)
for ii = 1:numIso
    ham = squeeze(ham_nuc(:,:,ii));
    [eigen_n, ~] = eig(ham);
    coef_n = randSphN(1, size(ham_nuc,1)); % Random configuration
    E_new(ii) = real(coef_n' * eigen_n' * ham * eigen_n * coef_n);

    Ispx = real(coef_n' * eigen_n' * Ix * eigen_n * coef_n)';
    Ispy = real(coef_n' * eigen_n' * Iy * eigen_n * coef_n)';
    Ispz = real(coef_n' * eigen_n' * Iz * eigen_n * coef_n)';
    newSpins(ii, :) = [Ispx, Ispy, Ispz];

    % Glauber algorithm to decide on updating the spin
    prob(ii) = 1 / (1 + exp((E_new(ii) - E_old(ii)) * beta));
end
selecSpins = prob >= rand("gpuArray"); % Random criteria for each isotope
nSpin(params.isoIdx(selecSpins), :) = newSpins(selecSpins, :);
