function coords = randSph(n, coordSystem)
    % Function to sample n points uniformly on a 3D spherical surface using Marsaglia method
    % Input: 
    %   n - Number of points to sample
    %   coordSystem - Coordinate system ('Cartesian' or 'Spherical')
    % Output:
    %   coords - n x 3 matrix containing coordinates of sampled points

    coords = zeros(n, 3);
    for i = 1:n
        u1 = 2*rand - 1; % uniformly distributed over [-1, 1]
        u2 = 2*rand - 1; % uniformly distributed over [-1, 1]
        s = u1^2 + u2^2;
        while s >= 1
            u1 = 2*rand - 1;
            u2 = 2*rand - 1;
            s = u1^2 + u2^2;
        end
        fac = 2*sqrt(1 - s);
        x = u1 * fac;
        y = u2 * fac;
        z = 1 - 2 * s;

        if strcmpi(coordSystem, 'cartesian')
            coords(i, :) = [x, y, z];
        elseif strcmpi(coordSystem, 'spherical')
            theta = atan2(y, x); % azimuthal angle
            phi = acos(z); % polar angle
            coords(i,:) = [theta, phi, 1]; % r = 1 for unit sphere
        else
            error('Invalid coordinate system specified. Use "Cartesian" or "Spherical".');
        end
    end
end
