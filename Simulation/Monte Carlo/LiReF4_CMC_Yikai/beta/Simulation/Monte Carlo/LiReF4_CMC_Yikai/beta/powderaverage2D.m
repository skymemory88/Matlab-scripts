function [FQav Qmod] = powderaverage2D(FQ)
%
% 2D powder average from simulated data
%
% Input: FQ is a 2D matrix with the calculated structure factors
%
% Output: FQav is a row vector containing the powder averaged structure
% factor and Qmod is the length of the scattering vector

% repeat in the plane to make powder average

Qx = -2:0.02:2;
Qx = repmat(Qx,length(Qx),1);

Qy = -2:0.02:2;
Qy = repmat(Qy',1,length(Qy));

FQ = repmat(FQ(1:end-1,1:end-1),4,4);
FQ(:,end+1) = FQ(:,1);
FQ(end+1,:) = FQ(1,:);  

% add third dimension

Qx = repmat(Qx,1,1,201);
Qy = repmat(Qy,1,1,201);
FQ = repmat(FQ,1,1,201);

z = -2:0.02:2;
for n = 1:length(z)
    Qz(:,:,n) = z(n)*ones(201,201);
end

% do the actual averaging
r = 0:0.035:2.035; % bins
for n = 1:length(r)-1  
    id{n} = (Qx.^2 + Qy.^2 + Qz.^2 >= r(n)^2) & (Qx.^2 + Qy.^2 + Qz.^2 < r(n+1)^2);
    FQav(n) = sum(FQ(id{n}))/(4*pi*(r(n)+0.035/2)^2);
end

Qmod = r(1:end-1)+0.035/2;

end

% Structural factor computation
function [kx ky FQ] = calc_mom_fft(r,m,Lcells)

[kx ky] = meshgrid(0:Lcells);
kz = 0.*kx;

SQx = zeros(size(kx));
SQy = zeros(size(kx));
SQz = zeros(size(kx));
FQ = zeros(size(kx));

for nk = 1:numel(kx)
    kk = [kx(nk) ky(nk) kz(nk)];
    for nr = 1:numel(r)/3
        SQx(nk) = SQx(nk) + exp(2i*pi*kk*r(:,nr))*m(1,nr);
        SQy(nk) = SQy(nk) + exp(2i*pi*kk*r(:,nr))*m(2,nr);
        SQz(nk) = SQz(nk) + exp(2i*pi*kk*r(:,nr))*m(3,nr);
    end
    
    Sperp = cross(kk/norm(kk),cross([SQx(nk) SQy(nk) SQz(nk)],kk/norm(kk)));
    FQ(nk) = dot(Sperp,Sperp);
end

FQ(isnan(FQ)) = 0;

end