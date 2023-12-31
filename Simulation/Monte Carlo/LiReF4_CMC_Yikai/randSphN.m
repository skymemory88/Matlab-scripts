function points = randSphN(N, dim)
% This function generates 'num_points' uniformly distributed points on an 'n'-dimensional hypersphere.

points = randn(N, dim);
norms = sqrt(sum(points.^2, 2));
points = points' ./ norms; % Normalizing each point to lie on the hypersphere
end