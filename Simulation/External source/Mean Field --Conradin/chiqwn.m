function Skw=chiqn(q,omega,epsilon,h,t)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.  

chi0=chi0n(h,t,omega,epsilon);

% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,size(q,1));

load('dipolwerte', 'Dpq', 'qvektor')

%q=qtemp;
for n=1:size(q,1)
   D(:,:,:,:,n)=dipole(q(n,:),12);
   dipolq(q(n,:));
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
  M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
  for n=1:N  %n,m Summe über Ionen in der Einheitszelle
    for m=1:N
      M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,nw)*D(:,:,n,m,nq);
    end
  end
  chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
    ((eye(size(M))-M)\([chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw);chi0(:,:,nw)]));
end
end



Skw=zeros(size(q,1),length(omega));

a=[5.175 0 0
   0 5.175 0
   0 0 10.75];
% Positions of the moments within the unit cell
tau=[0 0 0
     0 1/2 1/4
     1/2 1/2 1/2
     1/2 0 3/4];

% Convert tau to a
tau=tau*a;

% Unit cell volume
vol=sum(a(1,:).*cross(a(2,:),a(3,:)));
% Reciprocal lattice unit vectors
b=[2*pi*cross(a(2,:),a(3,:))/vol
   2*pi*cross(a(3,:),a(1,:))/vol
   2*pi*cross(a(1,:),a(2,:))/vol];
% Convert q from Miller indicies to reciprocal aangstroms


for nq=1:size(q,1)  
    % Convert q from Miller indicies to reciprocal aangstroms
    qb=q(nq,:)*b;
    for nw=1:length(omega)
        if t==0
            Skw(nq,nw)=1/pi* ...
                sum(sum((eye(3)-(qb'*qb)/(qb*qb')).* ...
                    imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
        else
            if omega(nw)~=0
                Skw(nq,nw)= 1/pi/(1-exp(-omega(nw)*11.6/t))* ...
                    sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                        imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
            else
                                Skw(nq,nw)=NaN;
            end

        end
    end
end


return
