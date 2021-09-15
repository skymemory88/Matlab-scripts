function chi=chi_qw(ion,q,chi0,dip_range,withdemagn,alpha)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.  


%% Initializes Dipole and Exchange(Zero for LiErF) at Each Q-point
N = size(ion.tau,1);
D_q=zeros(3,3,N,N,size(q,1));
E_q=zeros(3,3,N,N,size(q,1));
for nq=1:size(q,1)
     [D_q(:,:,:,:,nq),E_q(:,:,:,:,nq)]= calc_dip_ex(ion,q(nq,:),dip_range,withdemagn,alpha); 
end

chi=zeros(3,3,size(chi0,4),size(q,1));

for nw=1:size(chi0,4) 
    for nq=1:size(q,1)
      M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle  
      for n=1:N  %n,m Summe �ber Ionen in der Einheitszelle
        for m=1:N
          M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0(:,:,n,nw)*(D_q(:,:,n,m,nq)+E_q(:,:,n,m,nq)); % potentially .* instead of * --2021.04.12
        end
      end
      chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
        ((eye(size(M))-M)\([chi0(:,:,1,nw);chi0(:,:,2,nw);chi0(:,:,3,nw);chi0(:,:,4,nw)]));
    end
end
                
return
