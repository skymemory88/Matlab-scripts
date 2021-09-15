function chi=chi_qw2(q,chi0,sample)

if round(sample.Erbium*1e8) == 1e8
    Jex = sample.exEr;
    gL = sample.gLande_Er;
end
if round(sample.Holmium*1e8) == 1e8
    Jex = sample.exHo;
    gL = sample.gLande_Ho;
end

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,size(q,1));

for n=1:size(q,1)
    D(:,:,:,:,n) = gL^2*0.05368*dipole_direct(q(n,:),sample.dip_range) + exchange(q(n,:),Jex,sample.ex_range);
end

chi=zeros(3,3,size(chi0,4),size(q,1));
for nw=1:size(chi0,4)
    x0 = squeeze(mean(chi0(:,:,:,nw),3));
    for nq=1:size(q,1)
        M = squeeze(sum(sum(D(:,:,:,:,nq),4),3)/4); % average over the four sites in the unit cell and avoid double counting
        M = x0.*M;
        chi(:,:,nw,nq) = (ones(size(M))-M).\x0;
    end
end

return
