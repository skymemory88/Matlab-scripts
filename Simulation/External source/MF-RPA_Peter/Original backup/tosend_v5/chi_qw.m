function chi=chi_qw(q,chi0,sample)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.

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

% load('dipolwerte', 'Dpq', 'qvektor') % load file dipolwerte.mat, which contain variables Dpq and qvektor

%q=qtemp;
for n=1:size(q,1)
    D(:,:,:,:,n)=dipole_direct(q(n,:),100)+exchange(q(n,:),Jex);
end

%clear 'Dpq' 'qvektor'

% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*gL^2*0.05368;%*vol;

chi=zeros(3,3,size(chi0,4),size(q,1));
for nw=1:size(chi0,4)
    for nq=1:size(q,1)
        M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
        for n=1:N  %n,m Summe ueber Ionen in der Einheitszelle
            for m=1:N
                M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,n,nw)*D(:,:,n,m,nq);
            end
        end
        chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
            ((eye(size(M))-M)\([chi0(:,:,1,nw);chi0(:,:,2,nw);chi0(:,:,3,nw);chi0(:,:,4,nw)]));
    end
end


return
