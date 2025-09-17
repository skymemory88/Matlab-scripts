function points = randSphN(N, dim)
% This function generates N uniformly distributed points on an 'dim'-dimensional hypersphere.
points = randn(N, dim);
norms = sqrt(sum(points.^2, 2));
points = points' ./ norms; % Normalizing each point to lie on the hypersphere
end