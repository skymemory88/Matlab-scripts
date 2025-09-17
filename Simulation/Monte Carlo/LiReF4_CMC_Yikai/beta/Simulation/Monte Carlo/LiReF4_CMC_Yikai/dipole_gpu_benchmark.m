function dipole_gpu_benchmark()
% Benchmark and compare dipSum_a (CPU) and dipSum_gpu_a (GPU)
close all;
clearvars;

%% Setup test parameters
test_sizes = [5 15 35];  % Different lattice sizes
n_trials = 5;  % Number of trials for timing
n_random_tests = 10;  % Number of random configurations to test

% Initialize results storage
results = struct();

%% Run tests for different system sizes
for size_idx = 1:length(test_sizes)
    N = test_sizes(size_idx);
    fprintf('\n=== Testing %dx%dx%d system ===\n', N, N, N);

    % Create test system using your lattice function
    [params, const, gpuState] = create_realistic_test_system(N);

    % Report system info
    fprintf('Total spins: %d\n', size(params.pos, 1));
    n_neighbors = cellfun(@length, params.nList);
    fprintf('Avg neighbors: %.1f (min: %d, max: %d)\n', ...
        mean(n_neighbors), min(n_neighbors), max(n_neighbors));

    % Test dipSum implementations
    fprintf('\nTesting dipSum implementations...\n');
    [dipsum_results] = test_dipSum_all(params, const, gpuState, n_random_tests, n_trials);

    % Store results
    results(size_idx).N = N;
    results(size_idx).n_spins = size(params.pos, 1);
    results(size_idx).avg_neighbors = mean(n_neighbors);
    results(size_idx).dipsum = dipsum_results;
end

%% Generate comparison plots
generate_plots(results, test_sizes);

end

function [params, const, gpuState] = create_realistic_test_system(N)
% Create a test system using the actual lattice generation

% Setup ion parameters (using Ho as example)
ion.abc = {[5.162 0 0; 0 5.162 0; 0 0 10.70]};  % Er lattice parameters
ion.tau = [0     0     0
           0    1/2   1/4
          1/2   1/2   1/2
          1/2    0    3/4];
ion.idx = 1;

% Setup parameters
params.dims = [N, N, N];
params.coord = 'spherical';  % Use spherical cutoff
params.domRng = max(params.dims);  % Domain cutoff

% Dipole cutoff: 7 unit cells as in your code
params.dpRng = 15 * max(diag(ion.abc{ion.idx}));  % ~75 Angstrom for Ho

% Generate lattice
params.pos = lattice(ion.abc{ion.idx}, ion.tau, params);

fprintf('  Generated lattice with %d positions\n', size(params.pos, 1));
fprintf('  Dipole cutoff: %.1f Angstrom\n', params.dpRng);

% Build neighbor lists
[params.nList, params.nVecs, params.nDists] = neighborList_a(params.pos, params.dpRng);

% Constants (from your code)
const.mu0 = 4e-7 * pi;
const.muB = 9.274e-24;
const.J2meV = 6.24151e+21;
const.kB = 1.3806e-23;
gLande = 1.2;  % Approximate for Ho
const.gfac = const.mu0/4/pi * (gLande * const.muB)^2 * const.J2meV;

% Initialize GPU state
params.temp = 0.1;
params.field = [0; 0; 0];
gpuState = GPUSpinState(params, 1, 1);
gpuState.initialize(params, const.gfac);
end

function results = test_dipSum_all(params, const, gpuState, n_tests, n_trials)
% Test CPU and GPU dipSum implementations

N = size(params.pos, 1);

% Initialize result arrays
errors_gpu = zeros(n_tests, 1);

% Test accuracy with random configurations
fprintf('Testing accuracy with %d random configurations...\n', n_tests);
for test = 1:n_tests
    % Random spin configuration
    eSpin = randn(N, 3);
    eSpin = eSpin ./ vecnorm(eSpin, 2, 2);  % Normalize

    % Update GPU state
    gpuState.updateSpins(eSpin, 1, 1);

    % CPU calculation (reference)
    E_cpu = dipSum_a(params, const.gfac, eSpin);

    % GPU calculation
    E_gpu = dipSum_gpu_a(gpuState, 1, 1);

    % Store relative errors
    if abs(E_cpu) > 1e-10
        errors_gpu(test) = abs(E_gpu - E_cpu) / abs(E_cpu);
    else
        errors_gpu(test) = abs(E_gpu - E_cpu);
    end
end

% Timing tests with representative configuration
fprintf('Running timing tests...\n');
eSpin = randn(N, 3);
eSpin = eSpin ./ vecnorm(eSpin, 2, 2);
gpuState.updateSpins(eSpin, 1, 1);

% Pre-build pair data for GPU version (one-time cost)
fprintf('Building GPU pair data structure...\n');
tic;
E_gpu_first = dipSum_gpu_a(gpuState, 1, 1);  % First call builds pair data
t_build = toc;
fprintf('Pair data built in %.2f seconds\n', t_build);

% CPU timing
cpu_times = zeros(n_trials, 1);
for trial = 1:n_trials
    tic;
    E_cpu = dipSum_a(params, const.gfac, eSpin);
    cpu_times(trial) = toc;
end

% GPU timing (with pre-built pair data)
wait(gpuDevice);  % Ensure GPU is ready
gpu_times = zeros(n_trials, 1);
for trial = 1:n_trials
    tic;
    E_gpu = dipSum_gpu_a(gpuState, 1, 1);
    wait(gpuDevice);
    gpu_times(trial) = toc;
end

% Store results
results.errors_gpu = errors_gpu;
results.max_error_gpu = max(errors_gpu);
results.mean_error_gpu = mean(errors_gpu);
results.median_error_gpu = median(errors_gpu);

results.cpu_time_mean = mean(cpu_times);
results.cpu_time_std = std(cpu_times);
results.gpu_time_mean = mean(gpu_times);
results.gpu_time_std = std(gpu_times);

results.speedup_gpu = mean(cpu_times) / mean(gpu_times);
results.build_time = t_build;

% Verify both methods give same result
fprintf('Energy values: CPU=%.6e, GPU=%.6e, Relative diff=%.2e\n', ...
    E_cpu, E_gpu, abs(E_gpu - E_cpu)/abs(E_cpu));
end

function generate_plots(results, test_sizes)
% Generate comparison plots

% Figure 1: Error Analysis
figure('Position', [100, 100, 1200, 500]);
hold on;
for i = 1:length(results)
    resDip(i) = results(i).dipsum;
    if resDip(i).max_error_gpu > 0
        errors = resDip(i).errors_gpu;
        errors(errors == 0) = 1e-16;  % Replace zeros for log plot
        histogram(log10(errors), 20, 'DisplayName', ...
            sprintf('%d³ (%d spins)', test_sizes(i), results(i).n_spins));
    end
end
xlabel('log_{10}(Relative Error)');
ylabel('Frequency');
title('dipSum\_a: GPU vs CPU Relative Errors');
legend('Location', 'best');
grid on;
xlim([-16, 0]);

% Figure 2: Performance Comparison
figure('Position', [100, 650, 1200, 500]);
% Extract data
N_spins = [results.n_spins];
dipsum_speedups = zeros(size(results));
dipsum_cpu_times = zeros(size(results));
dipsum_gpu_times = zeros(size(results));
for i = 1:length(results)
    dipsum_speedups(i) = [resDip(i).speedup_gpu];
    dipsum_cpu_times(i) = [resDip(i).cpu_time_mean] * 1000;
    dipsum_gpu_times(i) = [resDip(i).gpu_time_mean] * 1000;
end

% Speedup plot
subplot(1,2,1);
hold on
loglog(N_spins, dipsum_speedups, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'dipSum\_a');
xlabel('Number of Spins');
ylabel('Speedup (CPU time / GPU time)');
title('GPU Speedup vs System Size');
legend('Location', 'best');
grid on;

% Add reference lines
loglog(N_spins, ones(size(N_spins)), 'k--', 'DisplayName', 'No speedup');

% Execution time plot
subplot(1,2,2);
hold on
loglog(N_spins, dipsum_cpu_times, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'dipSum CPU');
loglog(N_spins, dipsum_gpu_times, 's--', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'dipSum GPU');
xlabel('Number of Spins');
ylabel('Execution Time (ms)');
title('Execution Time vs System Size');
legend('Location', 'best');
grid on;

% Figure 3: Error Summary
figure('Position', [100, 200, 800, 400]);
hold on
% Extract errors
dipsum_max_errors = [resDip.max_error_gpu];
dipsum_mean_errors = [resDip.mean_error_gpu];

% Bar plot
x = 1:length(test_sizes);
width = 0.35;

semilogy(x + width/2, dipsum_max_errors, 'o', 'MarkerSize', 8, 'DisplayName', 'dipSum max');
semilogy(x + width/2, dipsum_mean_errors, 's', 'MarkerSize', 8, 'DisplayName', 'dipSum mean');

set(gca, 'XTick', x);
set(gca, 'XTickLabel', arrayfun(@(n) sprintf('%d³', n), test_sizes, 'UniformOutput', false));
xlabel('Lattice Size');
ylabel('Relative Error');
title('Maximum and Mean Relative Errors (CPU vs GPU)');
legend('Location', 'best');
grid on;
ylim([1e-16, 1e-6]);

%% Print summary statistics
fprintf('\n=== Accuracy Summary ===\n');
for i = 1:length(results)
    fprintf('\n%d³ lattice (%d spins, %.0f avg neighbors):\n', ...
        test_sizes(i), results(i).n_spins, results(i).avg_neighbors);
    fprintf('  dipSum_a:  max error = %.2e, mean = %.2e, median = %.2e\n', ...
        resDip(i).max_error_gpu, results(i).dipsum.mean_error_gpu, ...
        resDip(i).median_error_gpu);
end

fprintf('\n=== Performance Summary ===\n');
for i = 1:length(results)
    fprintf('\n%d³ lattice (%d spins):\n', test_sizes(i), results(i).n_spins);
    fprintf('  dipSum_a:  %.2fx speedup (CPU: %.3f ms, GPU: %.3f ms)\n', ...
        resDip(i).speedup_gpu, dipsum_cpu_times(i), dipsum_gpu_times(i));  % Fixed indexing
end

% Performance projection for 21x21x21
if test_sizes(end) < 21
    fprintf('\n=== Projected Performance for 21³ System ===\n');
    % Use power law fit for projection
    p_dipsum = polyfit(log(N_spins), log(dipsum_speedups), 1);

    N_21 = 21^3 * 4;  % 4 atoms per unit cell
    speedup_dipsum_21 = exp(polyval(p_dipsum, log(N_21)));

    fprintf('Projected speedups for 21³ system (%d spins):\n', N_21);
    fprintf('  dipSum_a:  %.1fx\n', speedup_dipsum_21);
end
end