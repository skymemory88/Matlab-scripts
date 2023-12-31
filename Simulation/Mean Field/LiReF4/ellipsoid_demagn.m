function N = ellipsoid_demagn(alpha)
% gets the demagnetization factor for a rotational ellipsis where a = b.
% alpha = c/a
% formulas are from 'Magnetic domains', A. Hubert and R, Schaefer, p. 121
if alpha == 1 %sphere
    Na = 1/3;
    Nb = 1/3;
    Nc = 1/3;
elseif alpha == 0 % needle
    Na = 0.5;
    Nb = 0.5;
    Nc = 0.0;
elseif alpha == inf % disk
    Na = 0;
    Nb = 0;
    Nc = 1;
elseif alpha < 1 % cigar
    Nc = (alpha^2/(1-alpha^2))*((1/sqrt(1-alpha^2))*asinh(sqrt(1-alpha^2)/alpha)-1);
    Nb = Nc;
    Na = 0.5*(1-Nc);
else % smarties
    Nc = (alpha^2/(alpha^2-1))*(1-(1/sqrt(alpha^2-1))*asin(sqrt(alpha^2-1)/alpha));
    Nb = Nc;
    Na = 0.5*(1-Nc);
end
N = [Na  0   0
      0  Nb  0
      0  0  Nc];
end