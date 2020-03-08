function [ion,evolution,E]=remf_pb(ion,h,t,withdemagn,alpha)

for i=1:size(ion.name,1)
    if ion.prop(i)~=0
        k=i;
    end
end

persistent dipole;
global rundipole

if(isempty(dipole))
    dipole=dipole_direct([0 0 0],80,ion.a{k});
    for i=1:size(ion.name,1)
        ion.cf(:,:,i)={cf(ion.J(i),ion.B(i,:))};
        ion.exch(:,:,i)={exchange([0,0,0],ion.ex(i))};
    end
end

if rundipole == true
    dipole=dipole_direct([0 0 0],80,ion.a{k});
    for i=1:size(ion.name,1)
        ion.cf(:,:,i)={cf(ion.J(i),ion.B(i,:))};
        ion.exch(:,:,i)={exchange([0,0,0],ion.ex(i))};
    end
    
    rundipole = false;
end

hvec=h;
momente_mean=0;
ion.mom_hyp=zeros([4,3,size(ion.name,1)]);
%%%
for i=1:size(ion.name,1)
    ion.mom_hyp(:,:,i)=ion.mom(:,:,i); %we will apply the virtual crystal approximation for the isotops of LiErF4
    momente_mean=momente_mean+ion.prop(i)*ion.mom(:,:,i);
end


%Lorenzterm = N/V * mu_0/4pi * (mu_b)^2 * 4pi/3
% Lorenz=3.124032/1000; %Ho
Lorenz=3.15452/1000; %Er

%gLande
for i=1:size(ion.name,1)
    ion.gLande(i)=gLande(ion.L(i),ion.S(i));
end

eins=zeros(3,3,4,4); eins(1,1,:,:)=1; eins(2,2,:,:)=1; eins(3,3,:,:)=1;
demagn=zeros(3,3,4,4);
demagn_t=ellipsoid_demagn(alpha);
demagn(1,1,:,:)=demagn_t(1,1);
demagn(2,2,:,:)=demagn_t(2,2);
demagn(3,3,:,:)=demagn_t(3,3);

%0.05368=mu_0*mu_B^2/4pi ([meV*Ang^-3], dipole_direct gives [Ang^3])
%0.000745836 = N/V * mu_0/4pi * (mu_b)^2

ion.d=zeros([3,3,4,4,size(ion.name,1)]);
ion.d_ex=zeros([3,3,4,4,size(ion.name,1)]);
for i=1:size(ion.name,1)
    ion.d(:,:,:,:,i)=ion.gLande(i)*(0.05368*dipole+eins*Lorenz/4-withdemagn*0.000745836*pi*demagn);
    ion.d_ex(:,:,:,:,i)=ion.exch{:,:,i};
    evolution.(ion.name{i})(1,:,:)=ion.mom(:,:,i);
    evolution.(ion.name_hyp{i})(1,:,:)=ion.mom_hyp(:,:,i);
end
alt=[1  1  1
    -1  1  1
    -1 -1  1
     1 -1  1];


% set convergence criteria. Convergence is accepted when either:
% (Difference < ConvMin) or (Difference < ConvMax and Niter>NiterMin) or (Niter>NiterMax)
% This ensures very good convergence when no problem, and decent
% convergence when difficult.
ConvMax=1e-6;
ConvMin=1e-7;
NiterMin=1000;
NiterMax=10000;

ion.mom_old=zeros([4,3,size(ion.name,1)]);
ion.mom_old_hyp=zeros([4,3,size(ion.name,1)]);

for iterations=1:NiterMax
 % iterations 
 for i=1:size(ion.name,1)
     ion.mom_old(:,:,i)=ion.mom(:,:,i);
     ion.mom_old_hyp(:,:,i)=ion.mom_hyp(:,:,i);
 end
 momente_old=momente_mean;
 E=zeros(1,4);
 for ionn=1:4
     %calculate meanfield
     h_dipol=zeros([size(ion.name,1),3]);
     h_dipol_hyp=zeros([size(ion.name,1),3]);
     h_ex=zeros([size(ion.name,1),3]);
     h_ex_hyp=zeros([size(ion.name,1),3]);
     for ionm=1:4
         for i=1:size(ion.name,1)
             h_dipol(i,:)=h_dipol(i,:)+ion.mom_old(ionm,:,i)*diag(ion.renorm(i,:))*ion.d(:,:,ionm,ionn,i)';
             h_ex(i,:)=h_ex(i,:)+ion.mom_old(ionm,:,i)*diag(ion.renorm(i,:))*ion.d_ex(:,:,ionm,ionn,i)';
             h_dipol_hyp(i,:)=h_dipol_hyp(i,:)+ion.mom_old_hyp(ionm,:,i)*diag(ion.renorm(i,:))*ion.d(:,:,ionm,ionn,i)';
             h_ex_hyp(i,:)=h_ex_hyp(i,:)+ion.mom_old_hyp(ionm,:,i)*diag(ion.renorm(i,:))*ion.d_ex(:,:,ionm,ionn,i)';
         end
     end
     
     %Virtual meanfield
     h_mf=zeros([size(ion.name,1),3]);
     for i=1:size(ion.name,1)
         h_mf(i,:)=ion.prop(i)*((1-ion.hyp(i))*(ion.gLande(i)*h_dipol(i,:)+h_ex(i,:)) + ion.hyp(i)*(ion.gLande(i)*h_dipol_hyp(i,:)+h_ex_hyp(i,:)));
         for j=1:size(ion.name,1)-1
             k=(i+j>size(ion.name,1))*size(ion.name,1);
             h_mf(i,:)=h_mf(i,:)+ion.prop(i+j-k)*((1-ion.hyp(i+j-k))*ion.gLande(i)*h_dipol(i+j-k,:)+ion.hyp(i+j-k)*ion.gLande(i)*h_dipol_hyp(i+j-k,:));%other kinds of ions
         end
     end

     %calculate moments of ions in a meanfield (ie diagonnalize the
     %hamiltonian)
%      energies=zeros(4,size(ion.name,1));
     for i=1:size(ion.name,1)
         if ion.prop(i)>0
             if ion.hyp(i)
%                  junk=[];
%                  [jx,jy,jz,energies(ionn,:)]=MF_moments_hyper(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i},ion.A(i),ion.I(i));
                 [jx,jy,jz,ix,iy,iz,energies(ionn,:)]=MF_moments_hyper(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i},ion.A(i),ion.I(i),ion.nuc(i));
                 ion.hyp_mom(ionn,:,i) = [ix iy iz];       % PB HACK
                 ion.mom_hyp(ionn,:,i)=update_moments([jx,jy,jz],evolution.(ion.name_hyp{i}),ionn,iterations);
             else
%                  junk=[];
                 [jx,jy,jz,energies(ionn,:)]=MF_moments(hvec,h_mf(i,:),t,ion.J(i),ion.gLande(i),ion.cf{:,:,i});
                 ion.hyp_mom(ionn,:,i) = [0 0 0];       % PB HACK
                 ion.mom(ionn,:,i)=update_moments([jx,jy,jz],evolution.(ion.name{i}),ionn,iterations);
             end
%               energies(ionn,:)=junk;
         end
     end
     E(ionn)=mean(energies(ionn,:));
 end
 momente_mean=0;
for i=1:size(ion.name,1)
    momente_mean=momente_mean+ion.prop(i)*((1-ion.hyp(i))*ion.mom(:,:,i)+ion.hyp(i)*ion.mom_hyp(:,:,i));
    evolution.(ion.name{i})(size(evolution.(ion.name{i}),1)+1,:,:)=ion.mom(:,:,i);
    evolution.(ion.name_hyp{i})(size(evolution.(ion.name_hyp{i}),1)+1,:,:)=ion.mom_hyp(:,:,i);
end

  
 % Iterate untill convergence criterion reached
  Diff=sum(sum(abs(momente_old-momente_mean)));
  if (Diff<ConvMin) || (Diff<ConvMax && iterations>NiterMin) 
    disp(num2str([t,h,iterations,mean(momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), n. %3.3i, <J>=(%3.3f,%3.3f,%3.3f)'))
    evolution.iterations=iterations;
    return
  end

end
disp(num2str([t,h,iterations,mean(momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), n. %3.3i, <J>=(%3.3f,%3.3f,%3.3f)'))
evolution.iterations=iterations;
return


function J=update_moments(Jnew,evolution,ionn,iterations)
global strategies;
if strategies.accelerator~=0 && iterations>10 && mod(iterations,strategies.expfit_period)~=0 % use accellerated update to shorten convergence
    Jold=squeeze(evolution(end-1,ionn,:))';
    Jnow=squeeze(evolution(end,ionn,:))';
    dnow=Jnow-Jold;
    J=Jnew+strategies.accelerator*dnow;
elseif strategies.expfit && mod(iterations,strategies.expfit_period)==0
    ndif=strategies.expfit_deltaN;
    Jold=squeeze(evolution(end-2*ndif+1,ionn,:))';
    Jnow=squeeze(evolution(end-1*ndif+1,ionn,:))';
    enNr=min(max(0,(Jnew-Jnow)./(Jnow-Jold)),0.998);
    if isempty(find(enNr>=1,1))
        feNNr=(Jnew-Jnow)./(enNr-1);
        J=Jnow-feNNr;
    else
        J=Jnew;
    end
else
    J=strategies.damping*squeeze(evolution(end,ionn,:))'+(1-strategies.damping)*Jnew;
end
return


function [jx,jy,jz,energies]=MF_moments(hvec,h_dipol,t,J,gLande,Hcf)
%Initiate J operators
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

%Calculate Hamiltonian
Hzeeman=(-gLande*0.05788)*(hvec(1)*Jx+hvec(2)*Jy+hvec(3)*Jz);
Hdipole=h_dipol(1)*Jx+h_dipol(2)*Jy+h_dipol(3)*Jz;
Ham=Hcf+Hzeeman-Hdipole;

%Diagonalize
[v,e]=eig(Ham);
E=e;
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);

%Calculate Matrixelements and Moments
if t==0
    % At zero temperature, use only lowest eigenvalue.
    energies=E(1,1);
    jx=real(v(:,1)'*Jx*v(:,1));
    jy=real(v(:,1)'*Jy*v(:,1));
    jz=real(v(:,1)'*Jz*v(:,1));
else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));
    jx=real(diag(v'*Jx*v))'*z;
    jy=real(diag(v'*Jy*v))'*z;
    jz=real(diag(v'*Jz*v))'*z;
    energies=sum(E)*z;
end

return

function [jx,jy,jz,ix,iy,iz,energies]=MF_moments_hyper(hvec,h_dipol,t,J,gLande,Hcf,A,I,nuc_mom)
% With hyperfine coupling
%Initiate J operators
Jz=diag(J:-1:-J);
Jzh=kron(Jz,eye(2*I+1));
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jph=kron(Jp,eye(2*I+1));
Jmh=kron(Jm,eye(2*I+1));
Jxh=(Jph+Jmh)/2;
Jyh=(Jph-Jmh)/2i;
%tensor product of cristal field to include nuclear moments
Hcfh=kron(Hcf,eye(2*I+1));
%Initiate I operators
Iz=diag(I:-1:-I);
Izh=kron(eye(2*J+1),Iz);
Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
Im=Ip';
Iph=kron(eye(2*J+1),Ip);
Imh=kron(eye(2*J+1),Im);
Ixh=(Iph+Imh)/2;
Iyh=(Iph-Imh)/2i;

%Calculate Hamiltonian
Hzeeman=(-gLande*0.05788)*(hvec(1)*Jxh+hvec(2)*Jyh+hvec(3)*Jzh);
Hzeeman_nuc=  (-3.152e-5*nuc_mom/I)*(hvec(1)*Ixh+hvec(2)*Iyh+hvec(3)*Izh);      %PB Hack

% Ham=Hcfh+Hzeeman-h_dipol(1)*Jxh-h_dipol(2)*Jyh-h_dipol(3)*Jzh + A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh);
Ham=Hcfh+Hzeeman+Hzeeman_nuc-h_dipol(1)*Jxh-h_dipol(2)*Jyh-h_dipol(3)*Jzh + A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh);

%Diagonalize
[v,e]=eig(Ham);
E=e;
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);

%Calculate Matrixelements and Moments
if t==0
    % At zero temperature, use only lowest eigenvalue.
    energies=E(1,1);
    jx=real(v(:,1)'*Jxh*v(:,1));
    jy=real(v(:,1)'*Jyh*v(:,1));
    jz=real(v(:,1)'*Jzh*v(:,1));
    % PB HACK
    ix=real(v(:,1)'*Ixh*v(:,1));
    iy=real(v(:,1)'*Iyh*v(:,1));
    iz=real(v(:,1)'*Izh*v(:,1));
    % END
else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));
    jx=real(diag(v'*Jxh*v))'*z;
    jy=real(diag(v'*Jyh*v))'*z;
    jz=real(diag(v'*Jzh*v))'*z;
    % PB HACK
    ix=real(diag(v'*Ixh*v))'*z;
    iy=real(diag(v'*Iyh*v))'*z;
    iz=real(diag(v'*Izh*v))'*z;
    % END
    energies=sum(E)*z;
end

return
