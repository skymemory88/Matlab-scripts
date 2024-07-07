function BraggPeaks(pos, spin, wavelength, kVec, Xaxis, xrange)
% Function to simulate neutron powder diffraction or single crystal diffraction data including magnetic peaks
% Parameters:
% pos: coordinate (angstrom) of the N ions in the crystal (N x 3)
% spin: N x 3 matrix of spin vectors (magnetic moments) corresponding to ion positions
% wavelength: Neutron wavelength in Angstroms
% kVec: Incident orientation specified as [kx, ky, kz] for single crystal diffraction
% Xaxis: 'space' for d-spacing, 'angle' for 2-theta
% xrange: Range for X-axis values

% Validate inputs
spin = squeeze(spin);
if size(pos, 2) ~= 3
    error('pos must be N x 3 matrices');
end
if ~isempty(kVec) && length(kVec) ~= 3
    error('kVec must be a vector of length 3');
end
if ~ismember(Xaxis, {'space', 'angle'})
    error('Xaxis must be either ''space'' or ''angle''');
end
if strcmp(Xaxis, 'angle')
    % Neutron Powder Diffraction
    if isempty(xrange)
        theta = linspace(0, 180, 361);
    else
        theta = linspace(xrange(1), xrange(2), 4*180);
    end
    pattern = calculate_diffraction(pos, spin, theta, wavelength, 'angle');

    % Plot the simulated diffraction pattern as a function of 2-theta
    figure;
    plot(theta, pattern, 'LineWidth', 1.5);
    xlabel('2$\theta$ (degrees)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Intensity (a.u.)', 'FontSize', 14);
    title('Simulated Neutron Powder Diffraction Pattern');
    grid on;
else
    % Single Crystal Neutron Diffraction
    pairDist = pdist2(pos, pos);
    max_extent = max(pairDist(:));
    q_min = 2 * pi / max_extent; % minimal q resolution
    % Calculate max_extent using pdist2 for the true maximum distance
    q_max = 4 * pi / wavelength; % Maximum q value based on the smallest d-spacing (wavelength / 2)
    if ~isempty(xrange)
        xrange(xrange == 0) = 0.01;
        xrange = 2*pi ./ xrange;
        if min(xrange) > q_min; q_min = min(xrange); end
        if max(xrange) < q_max; q_max = max(xrange); end
    end
    q_values = linspace(q_min, q_max, 1000); % Range of q values
    pattern = calculate_diffraction(pos, spin, q_values, kVec, 'space');

    % Calculate d-spacing corresponding to q_values
    d_range = 2*pi ./ q_values;

    % Plot the simulated diffraction pattern as a function of d-spacing
    figure;
    plot(d_range, pattern, 'LineWidth', 1.5);
    xlabel('d-spacing (\AA)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Intensity (a.u.)', 'FontSize', 14);
    title(sprintf('Simulated Neutron Single Crystal Diffraction [%d %d %d]', kVec(1), kVec(2), kVec(3)));
    grid on;
end
end

function pattern = calculate_diffraction(pos, spin, values, direction, mode)
% Calculate the diffraction pattern for both powder and single crystal
% Parameters:
% pos: N x 3 matrix of ion positions (N ions in the crystal)
% spin: N x 3 matrix of spin vectors (magnetic moments)
% values: q_values or theta values depending on mode
% direction: kVec for single crystal, wavelength for powder
% mode: 'angle' for powder diffraction, 'space' for single crystal diffraction

num_values = length(values);
pattern = zeros(num_values, 1);
if strcmp(mode, 'angle')
    % For powder diffraction, use 2-theta values
    wavelength = direction; % direction stores wavelength in this case
    two_theta = values;
    q_values = 4 * pi * sind(two_theta / 2) / wavelength;
    q_vectors = q_values(:) .* [1, 1, 1]; % Ensure q_vector has the correct dimensions
    phase_factors = exp(2 * pi * 1i * (pos * q_vectors')); % pos * q_vector' to match dimensions
    for ii = 1:num_values
        q_vector = q_vectors(ii, :);
        F_ionic = sum(phase_factors(:, ii)); % Ionic structure factor
        if ~isempty(spin)
            % Calculate the magnetic structure factor
            F_magnetic = sum(phase_factors(:, ii) .* vecnorm(spin - sum(spin .* q_vector, 2) .* (q_vector/norm(q_vector)), 2, 2)); % Magnetic structure factor
        else
            F_magnetic = 0;
        end
        F_q = F_ionic + F_magnetic; % Combined structure factor
        pattern(ii) = abs(F_q)^2;
    end
else
    % For single crystal diffraction, use q_values and kVec
    q_values = values;
    kVec = direction; % propogation vector of the incoming neutron beam
    q_vectors = q_values(:) .* kVec; % Ensure correct dimensions for q_vectors
    phase_factors = exp(2 * pi * 1i * (pos * q_vectors')); % Calculate phase factors for all q-values at once
    for ii = 1:num_values
        q_vector = q_vectors(ii, :);
        F_ionic = sum(phase_factors(:, ii)); % Ionic structure factor
        if ~isempty(spin)
            % Calculate the magnetic structure factor
            F_magnetic = sum(phase_factors(:, ii) .* vecnorm(spin - sum(spin .* q_vector, 2) .* (q_vector/norm(q_vector)), 2, 2)); % Magnetic structure factor
        else
            F_magnetic = 0;
        end
        F_q = F_ionic + F_magnetic; % Combined structure factor
        pattern(ii) = abs(F_q)^2;
    end
end

% Normalize pattern
if ~isempty(pattern)
    pattern = pattern / max(pattern) * 100;
end
end
