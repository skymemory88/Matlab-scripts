function coords = randSph_GPU(n, coordSystem)
    % Function to sample n points uniformly on a 3D spherical surface using Marsaglia method
    % Input: 
    %   n - Number of points to sample
    %   coordSystem - Coordinate system ('Cartesian' or 'Spherical')
    % Output:
    %   coords - n x 3 matrix containing coordinates of sampled points

    % Pre-allocate GPU arrays
    coords = gpuArray.zeros(n, 3);
    valid = gpuArray.zeros(n, 1, 'logical');

    while ~all(valid)
        % Generate points within the unit square
        u1 = 2*gpuArray.rand(n, 1, 'single') - 1; % uniformly distributed over [-1, 1]
        u2 = 2*gpuArray.rand(n, 1, 'single') - 1; % uniformly distributed over [-1, 1]
        s = u1.^2 + u2.^2;

        % Identify points within the unit circle to ensure uniform distribution on sphere
        withinCircle = s < 1;
        newValid = withinCircle & ~valid;  % New valid points not previously validated

        % Calculate coordinates only for valid points
        fac = 2.*sqrt(1 - s(newValid));
        x = u1(newValid) .* fac;
        y = u2(newValid) .* fac;
        z = 1 - 2 * s(newValid);

        % Store the new valid coordinates
        coords(newValid, :) = [x, y, z];
        valid(newValid) = true;
    end

    % Convert coordinates from Cartesian to Spherical if necessary
    if strcmpi(coordSystem, 'spherical')
        x = coords(:, 1);
        y = coords(:, 2);
        z = coords(:, 3);
        theta = atan2(y, x); % azimuthal angle
        phi = acos(z); % polar angle
        coords = [theta, phi, ones(n, 1, 'like', theta)]; % r = 1 for unit sphere
    end

    % Optionally, gather the results back to CPU if subsequent operations require CPU processing
    % coords = gather(coords);
end
