function chi=chi_qw(q,omega,epsilon,h,t)

chi0=zeros(3,3,length(omega));
for n=1:length(omega)
  chi0(:,:,n)=chi0_w(h,t,omega(n),epsilon);
end

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,length(omega));
load dipole_saves
qz=qz;
for n=1:size(q,1)
  mz=find((q(n,1)==qz(:,1)) & (q(n,2)==qz(:,2)) & (q(n,3)==qz(:,3)));
  my=find((q(n,1)==qy(:,1)) & (q(n,2)==qy(:,2)) & (q(n,3)==qy(:,3)));
  mx=find((q(n,1)==qx(:,1)) & (q(n,2)==qx(:,2)) & (q(n,3)==qx(:,3)));
  if ~isempty(mx)
    D(:,:,:,:,n)=dqx(:,:,:,:,mx);
  elseif ~isempty(my)
    D(:,:,:,:,n)=dqy(:,:,:,:,my);
  elseif ~isempty(mz)
    D(:,:,:,:,n)=dqz(:,:,:,:,mz);
  else
    D(:,:,:,:,n)=dipole(q(n,:),12);
  end
end
% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*(5/4)^2*0.05368;%*vol;

chi=zeros(3,3,length(omega),size(q,1));
for nw=1:length(omega)
for nq=1:size(q,1)
  M=zeros(3*N);
  for n=1:N
    for m=1:N
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,nw)*D(:,:,n,m,nq);
    end
  end
  chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw)]));
end
end


return
