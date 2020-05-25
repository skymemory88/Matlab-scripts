function [alpha, beta]=random_angles
% gives two random angles alpha and beta to give in input to cl_eff_mod,
% with a uniform distribution of the vector state = cos(alpha)|+> +
% exp(1i*beta)sin(alpha)|-> on the Bloch sphere.

% uses Marsaglia's method to pick a random vector on the sphere (cf
% Marsaglia, G. (1972). "Choosing a Point from the Surface of a Sphere".
% Ann. Math. Stat. 43 (2): 645–646. )
x1=2*rand(1)-1;
x2=2*rand(1)-1;
while(x1^2+x2^2>=1)
    x1=2*rand(1)-1;
    x2=2*rand(1)-1;
end

Vbloch=[2*x1*sqrt(1-x1^2-x2^2) ; 2*x2*sqrt(1-x1^2-x2^2) ; 1-2*(x1^2+x2^2)];

% finds alpha and beta by computing the spherical coordinates of Vbloch

alpha=acos(Vbloch(3));

beta=acos(Vbloch(1)/sin(alpha));
s=asin(Vbloch(2)/sin(alpha));
beta=sign(s)*beta;

end