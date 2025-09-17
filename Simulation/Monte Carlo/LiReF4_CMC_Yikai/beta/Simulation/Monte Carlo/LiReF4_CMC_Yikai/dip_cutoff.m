function [optimal_cutoff, results] = dip_cutoff(ion, params_base, const)
% Analyze dipole-dipole interaction energy vs cutoff distance
% Tests different lattice sizes to find optimal cutoff in units of unit cells
%
% Input:
%   ion - ion structure with .abc (lattice parameters) and .tau (basis)
%   params_base - base parameters (will be modified with different dims)
%   const - constants including .gfac
%
% Output:
%   optimal_cutoff - recommended cutoff in unit cells
%   results - detailed analysis results

%% Setup
% Get lattice parameters
abc = ion.abc{ion.idx};
tau = ion.tau;
lattice_const = min(diag(abc)); % Minimum lattice constant

% Test different system sizes (in unit cells)
dims_range = 5:5:50; % From 3x3x3 to 15x15x15
n_tests = length(dims_range);

% Initialize storage
E_total_vs_size = zeros(n_tests, 1);
E_vs_distance = cell(n_tests, 1);
r_bins = cell(n_tests, 1);
E_accumulated = cell(n_tests, 1);
r_accumulated = cell(n_tests, 1);

%% Generate test spin configuration
% Random spins normalized to unit length
rng(42); % For reproducibility
generate_test_spins = @(N) bsxfun(@rdivide, randn(N, 3), vecnorm(randn(N, 3), 2, 2));

for idx = 1:n_tests
    dims = dims_range(idx) * [1 1 1];
    
    % Create lattice
    params = params_base;
    params.dims = dims;
    params.coord = 'cubic'; % Use cubic shape for clear cutoff analysis
    params.domRng = max(dims);
    
    tic;
    pos = lattice(abc, tau, params);
    N = size(pos, 1);
    
    % Generate spins
    % spin = generate_test_spins(N);
    spin = repmat([0 0 1], N, 1);
    
    % Find central site
    pos0 = mean(pos, 1);
    dist0 = vecnorm(pos - pos0, 2, 2);
    [~, idx0] = min(dist0);
    
    % Calculate distances from center
    r_from_center = vecnorm(pos - pos(idx0, :), 2, 2);
    r_from_center(idx0) = []; % Remove self
    
    % Calculate dipole energy contribution from each spin
    r_vec = pos - pos(idx0, :);
    r_vec(idx0, :) = [];
    r_vec = r_vec * 1e-10; % Convert to meters
    dist = vecnorm(r_vec, 2, 2);
    
    spin_center = spin(idx0, :);
    spin_others = spin;
    spin_others(idx0, :) = [];
    
    % Vectorized dipole calculation
    ss_dot = spin_others * spin_center';
    sr_dot1 = r_vec * spin_center' ./ dist;
    sr_dot2 = sum(r_vec .* spin_others, 2) ./ dist;
    
    E_individual = const.gfac * (ss_dot ./ dist.^3 - 3 * (sr_dot1 .* sr_dot2) ./ dist.^5);
    
    % Store total energy
    E_total_vs_size(idx) = sum(E_individual);
    
    % Bin by distance for analysis
    dr = lattice_const / 2; % Bin width
    r_max = max(r_from_center);
    r_edges = 0:dr:r_max+dr;
    r_centers = (r_edges(1:end-1) + r_edges(2:end)) / 2;
    
    % Calculate energy in each distance bin
    E_binned = zeros(length(r_centers), 1);
    for i = 1:length(r_centers)
        mask = (r_from_center >= r_edges(i)) & (r_from_center < r_edges(i+1));
        E_binned(i) = sum(E_individual(mask));
    end
    
    % Store results
    E_vs_distance{idx} = E_binned;
    r_bins{idx} = r_centers;
    
    % Calculate accumulated energy vs cutoff
    r_cutoff_test = linspace(0, r_max, 200);
    E_acc = zeros(size(r_cutoff_test));
    for i = 1:length(r_cutoff_test)
        mask = r_from_center <= r_cutoff_test(i);
        E_acc(i) = sum(E_individual(mask));
    end
    E_accumulated{idx} = E_acc;
    r_accumulated{idx} = r_cutoff_test;
    
    t_calc = toc;
    fprintf('%17d | %10d | %13.6f | %8.3f\n', dims(1), N, E_total_vs_size(idx), t_calc);
end

%% Analyze convergence
fprintf('\n\nAnalyzing convergence...\n');

% Find where energy converges
convergence_threshold = 0.99; % 99% of converged value
cutoff_unit_cells = zeros(n_tests, 1); % Convert to unit cells

for idx = 1:n_tests
    % Find cutoff where 99% of energy is captured
    E_target = convergence_threshold * E_accumulated{idx}(end);
    cutoff_idx = find(E_accumulated{idx} >= E_target, 1);
    if ~isempty(cutoff_idx)
        cutoff_unit_cells(idx) = r_accumulated{idx}(cutoff_idx) / lattice_const;
    end
end

% Determine optimal cutoff
% Use the point where cutoff stops increasing with system size
diff_cutoff = diff(cutoff_unit_cells);
stable_idx = find(abs(diff_cutoff) < 0.1, 1);
if isempty(stable_idx)
    stable_idx = length(cutoff_unit_cells) - 1;
end
optimal_cutoff = cutoff_unit_cells(stable_idx);

%% Create figures
% Figure 1: Energy vs distance (radial distribution)
figure('Name', 'Dipole Energy vs Distance', 'Position', [100 100 1200 500]);

subplot(1, 2, 1);
hold on;
colors = parula(n_tests);
for idx = 1:n_tests
    plot(r_bins{idx}/lattice_const, E_vs_distance{idx}, '-', ...
         'Color', colors(idx,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%dx%dx%d', dims_range(idx), dims_range(idx), dims_range(idx)));
end
xlabel('Distance (unit cells)');
ylabel('Energy contribution (meV)');
title('Dipole-Dipole Energy vs Distance');
legend('Location', 'eastoutside');
grid on;
set(gca, 'YScale', 'log');
xlim([0 max(dims_range)/2]);

subplot(1, 2, 2);
hold on;
for idx = 1:n_tests
    plot(r_bins{idx}/lattice_const, abs(E_vs_distance{idx}), '-', ...
         'Color', colors(idx,:), 'LineWidth', 1.5);
end
xlabel('Distance (unit cells)');
ylabel('|Energy contribution| (meV)');
title('Absolute Energy vs Distance (log scale)');
grid on;
set(gca, 'YScale', 'log');
xlim([0 max(dims_range)/2]);

% Figure 2: Accumulated energy vs cutoff
figure('Name', 'Accumulated Energy vs Cutoff', 'Position', [100 650 1200 500]);

subplot(1, 2, 1);
hold on;
for idx = 1:n_tests
    plot(r_accumulated{idx}/lattice_const, E_accumulated{idx}, '-', ...
         'Color', colors(idx,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%dx%dx%d', dims_range(idx), dims_range(idx), dims_range(idx)));
end
xlabel('Cutoff distance (unit cells)');
ylabel('Accumulated energy (meV)');
title('Total Energy within Cutoff Radius');
legend('Location', 'southeast');
grid on;
xlim([0 max(dims_range)/2]);

subplot(1, 2, 2);
hold on;
for idx = 1:n_tests
    E_fraction = E_accumulated{idx} / E_accumulated{idx}(end);
    plot(r_accumulated{idx}/lattice_const, E_fraction * 100, '-', ...
         'Color', colors(idx,:), 'LineWidth', 2);
end
xlabel('Cutoff distance (unit cells)');
ylabel('Fraction of total energy (%)');
title('Energy Convergence vs Cutoff');
grid on;
xlim([0 max(dims_range)/2]);
ylim([0 105]);

% Add reference lines
plot(xlim, [95 95], 'k--', 'LineWidth', 1);
plot(xlim, [99 99], 'k--', 'LineWidth', 1);

% Figure 3: Convergence analysis
figure('Name', 'Convergence Analysis', 'Position', [100 100 800 600]);

subplot(2, 2, 1);
plot(dims_range, E_total_vs_size, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('System size (unit cells per side)');
ylabel('Total energy (meV)');
title('Total Energy vs System Size');
grid on;

subplot(2, 2, 2);
plot(dims_range, cutoff_unit_cells, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('System size (unit cells per side)');
ylabel('99% energy cutoff (unit cells)');
title('Convergence of Cutoff Distance');
grid on;

subplot(2, 2, 3);
% Energy per shell analysis for largest system
shell_width = 0.5; % in unit cells
max_shells = floor(max(dims_range)/2 / shell_width);
shell_energy = zeros(max_shells, 1);
shell_count = zeros(max_shells, 1);

for i = 1:max_shells
    r_min = (i-1) * shell_width * lattice_const;
    r_max = i * shell_width * lattice_const;
    mask = (r_from_center >= r_min) & (r_from_center < r_max);
    shell_energy(i) = sum(abs(E_individual(mask)));
    shell_count(i) = sum(mask);
end

bar(1:max_shells, shell_energy);
xlabel('Shell number (0.5 unit cell width)');
ylabel('|Energy| per shell (meV)');
title('Energy Contribution by Distance Shell');
set(gca, 'YScale', 'log');
grid on;

subplot(2, 2, 4);
% Cumulative energy fraction
cumulative_fraction = cumsum(shell_energy) / sum(shell_energy);
plot(1:max_shells, cumulative_fraction * 100, 'g-', 'LineWidth', 2);
xlabel('Shell number (0.5 unit cell width)');
ylabel('Cumulative energy fraction (%)');
title('Cumulative Energy by Shell');
grid on;
ylim([0 105]);

%% Summary output
fprintf('\n=== ANALYSIS SUMMARY ===\n');
fprintf('Lattice constant: %.2f Å\n', lattice_const);
fprintf('Tested system sizes: %d×%d×%d to %d×%d×%d unit cells\n', ...
        dims_range(1), dims_range(1), dims_range(1), ...
        dims_range(end), dims_range(end), dims_range(end));

%% Store results
results.unit_cells_tested = dims_range;
results.E_total_vs_size = E_total_vs_size;
results.E_vs_distance = E_vs_distance;
results.r_bins = r_bins;
results.E_accumulated = E_accumulated;
results.r_accumulated = r_accumulated;
results.cutoff_unit_cells = cutoff_unit_cells;

end

%% Lattice generation function (simplified version)
function pos = lattice(uVec, bVec, params)
    % Construct a lattice with a finite size and custom cutoff shape
    % Input:
    %   bVec - n x 3 matrix containing the basis vectors as rows
    %   uVec - n x 3 matrix containing the unit vectors as rows
    %   params.dims - Vector [nx, ny, nz] specifying the number of unit cells
    %   params.domRng - domain range for cutoff
    %   params.coord - coordinate system type ('cubic' or 'spherical')
    
    dims = params.dims;
    cutR = params.domRng/2 * min(diag(uVec)); % lattice cutoff radius
    maskType = params.coord;
    
    % Generate unit cell indices centered around zero
    [X, Y, Z] = ndgrid(-floor(dims(1)/2):ceil(dims(1)/2)-1, ...
                       -floor(dims(2)/2):ceil(dims(2)/2)-1, ...
                       -floor(dims(3)/2):ceil(dims(3)/2)-1);
    
    % Reshape to vectors
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    % Generate base positions for all unit cells
    basePos = [X, Y, Z] * uVec;
    
    % Correct the unit vector as they are defined by ratio instead of absolute length
    bVec = bVec * uVec;
    
    % Generate all lattice points based on basis vectors
    allPos = [];
    for i = 1:size(bVec, 1)
        allPos = [allPos; basePos + repmat(bVec(i,:), size(basePos, 1), 1)];
    end
    
    % Apply cutoff
    if strcmp(maskType, 'spherical')
        cutoffMask = vecnorm(allPos, 2, 2) <= cutR;
    else % cubic
        cutoffMask = all(abs(allPos) <= cutR, 2);
    end
    
    % Final positions
    pos = allPos(cutoffMask, :);
end