function [momente_er,momente_er_hyp,momente_ho,evolution]=remf(h,t,momente_er,momente_ho,exHo,exEr,xEr,xHo,ErHyp,HoHyp,withdemagn,alpha,renorm_Er,renorm_Ho)

persistent dipole cryst_Ho cryst_Er exch_Er exch_Ho;

if(isempty(dipole))
    dipole=dipole_direct([0 0 0],100);
    cryst_Ho=cf_Ho;
    cryst_Er=cf_Er;
    exch_Er=exchange([0,0,0],exEr);
    exch_Ho=exchange([0,0,0],exHo);
end

hvec=h;

%%%

momente_er_hyp=momente_er; %we will apply the virtual crystal approximation for the isotops of LiErF4
momente_mean=xEr*momente_er+xHo*momente_ho;

%Lorenzterm = N/V * mu_0/4pi * (g_l*mu_b)^2 * 4pi/3
Lorenz=3.124032/1000;

%J and gLande
J_Er=15/2;
gLande_Er=1.2;
J_Ho=8;
gLande_Ho=1.25;

%interaction matrix, d_RE(:,:,i,j) is the 3by3 interaction matrix between
%site i and j.
eins=zeros(3,3,4,4); eins(1,1,:,:)=1; eins(2,2,:,:)=1; eins(3,3,:,:)=1;
demagn=zeros(3,3,4,4);
demagn_t=ellipsoid_demagn(alpha);
demagn(1,1,:,:)=demagn_t(1,1);
demagn(2,2,:,:)=demagn_t(2,2);
demagn(3,3,:,:)=demagn_t(3,3);

%0.05368=mu_0*mu_B^2/4pi ([meV*Ang^-3], dipole_direct gives [Ang^3])
%0.000745836 = N/V * mu_0/4pi
% d_er=(gLande_Er)^2*0.05368*dipole_direct([0 0 0],15)+exchange([0,0,0],Jex_er)+eins*Lorenz*(gLande_Er)^2/4-(0.000745836*(gLande_Er)^2)*demagn/4;
% d_ho=(gLande_Ho)^2*0.05368*dipole_direct([0 0 0],15)+exchange([0,0,0],Jex_ho)+eins*Lorenz*(gLande_Ho)^2/4-(0.000745836*(gLande_Ho)^2)*demagn/4;
d_er=gLande_Er*(0.05368*dipole+eins*Lorenz/4-withdemagn*0.000745836*pi*demagn);
d_ho=gLande_Ho*(0.05368*dipole+eins*Lorenz/4-withdemagn*0.000745836*pi*demagn);
d_ex_er=exch_Er;
d_ex_ho=exch_Ho;
alt=[1 1 1
    -1 -1 1
    -1 -1 1
    1 1 1];

evolution.ho(1,:,:)=momente_ho;
evolution.er(1,:,:)=momente_er;
evolution.er_hyp(1,:,:)=momente_er;

% set convergence criteria. Convergence is accepted when either:
% (Difference < ConvMin) or (Difference < ConvMax and Niter>NiterMin) or (Niter>NiterMax)
% This ensures very good convergence when no problem, and decent
% convergence when difficult.
ConvMax=1e-6;
ConvMin=1e-7;
NiterMin=1000;
NiterMax=10000;

for iterations=1:NiterMax
 % iterations 
 momente_old_er=momente_er;
 momente_old_er_hyp=momente_er_hyp;
 momente_old_ho=momente_ho;
 momente_old=momente_mean;
 
 for ionn=1:4
     %calculate meanfield
     h_dipol_er=[0 0 0];
     h_ex_er=[0 0 0];
     h_dipol_er_hyp=[0 0 0];
     h_ex_er_hyp=[0 0 0];
     h_dipol_ho=[0 0 0];
     h_ex_ho=[0 0 0];
     for ionm=1:4
         h_dipol_er=h_dipol_er+momente_old_er(ionm,:)*diag(renorm_Er)*d_er(:,:,ionm,ionn)';
         h_ex_er=h_ex_er+momente_old_er(ionm,:)*diag(renorm_Er)*d_ex_er(:,:,ionm,ionn)';
         h_dipol_er_hyp=h_dipol_er_hyp + momente_old_er_hyp(ionm,:)*diag(renorm_Er)*d_er(:,:,ionm,ionn)';
         h_ex_er_hyp=h_ex_er_hyp + momente_old_er_hyp(ionm,:)*diag(renorm_Er)*d_ex_er(:,:,ionm,ionn)';
         h_dipol_ho=h_dipol_ho + momente_old_ho(ionm,:)*diag(renorm_Ho)*d_ho(:,:,ionm,ionn)';
         h_ex_ho=h_ex_ho + momente_old_ho(ionm,:)*diag(renorm_Ho)*d_ex_ho(:,:,ionm,ionn)';
     end
     %h_mf is the virtual crystal mean field.
     h_mf_ho=xEr*((1-ErHyp)*gLande_Ho*h_dipol_er + ErHyp*gLande_Ho*h_dipol_er_hyp)+xHo*(gLande_Ho*h_dipol_ho+h_ex_ho);
     h_mf_er=xEr*((1-ErHyp)*(gLande_Er*h_dipol_er+h_ex_er) + ErHyp*(gLande_Er*h_dipol_er_hyp+h_ex_er_hyp))+xHo*gLande_Er*h_dipol_ho;

     %calculate moments of ions in a meanfield (ie diagonnalize the
     %hamiltonian)
     if xEr>0
         if ErHyp > 0
             junk=[];
             A_Er=0.00043412; %hyperfine coupling energy in meV for Erbium (Phys. Rev. B 2, 2298 - 2301 (1970))
             I_Er=3.5;
             [jx,jy,jz,junk]=MF_moments_hyper(hvec,h_mf_er,t,J_Er,gLande_Er,cryst_Er,A_Er,I_Er);
             momente_er_hyp(ionn,:)=update_moments([jx,jy,jz],evolution.er_hyp,ionn,iterations);
         end
         if ErHyp~=1
             junk=[];
             [jx,jy,jz,junk]=MF_moments(hvec,h_mf_er,t,J_Er,gLande_Er,cryst_Er);
             momente_er(ionn,:)=update_moments([jx,jy,jz],evolution.er,ionn,iterations);
         end
     end
     if xHo>0
         if HoHyp
           junk=[];
           A_Ho=0.003361; %hyperfine coupling energy in meV (from PRB 75 054426).
           I_Ho=3.5;
           [jx,jy,jz,junk]=MF_moments_hyper(hvec,h_mf_ho,t,J_Ho,gLande_Ho,cryst_Ho,A_Ho,I_Ho);
           momente_ho(ionn,:)=update_moments([jx,jy,jz],evolution.ho,ionn,iterations);
         else
           junk=[];
           [jx,jy,jz,junk]=MF_moments(hvec,h_mf_ho,t,J_Ho,gLande_Ho,cryst_Ho);
           momente_ho(ionn,:)=update_moments([jx,jy,jz],evolution.ho,ionn,iterations);
         end
     end
 end

  momente_mean=xEr*((1-ErHyp)*momente_er+ErHyp*momente_er_hyp)+xHo*momente_ho;
  %momente_mean=(xHo*momente_ho+xEr*momente_er)/(xHo+xEr);
  
  evolution.ho(size(evolution.ho,1)+1,:,:)=momente_ho;
  evolution.er(size(evolution.er,1)+1,:,:)=momente_er;
  evolution.er_hyp(size(evolution.er_hyp,1)+1,:,:)=momente_er_hyp;
  
 % Iterate untill convergence criterion reached
  Diff=sum(sum(abs(momente_old-momente_mean)));
  if (Diff<ConvMin) || (Diff<ConvMax && iterations>NiterMin) 
    disp(num2str([t,h,iterations,mean(alt.*momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f)'))
    evolution.iterations=iterations;
    return
  end

end
disp(num2str([t,h,iterations,mean(alt.*momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f)'))
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


function [jx,jy,jz,e]=MF_moments(hvec,h_dipol,t,J,gLande,Hcf)
%Initiate J operators
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

%Calculate Hamiltonian
Hzeeman=(-gLande*0.05788)*(hvec(1)*Jx+hvec(2)*Jy+hvec(3)*Jz);
Ham=Hcf+Hzeeman-h_dipol(1)*Jx-h_dipol(2)*Jy-h_dipol(3)*Jz;

%Diagonalize
[v,e]=eig(Ham);
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);

%Calculate Matrixelements and Moments
if t==0
    % At zero temperature, use only lowest eigenvalue.
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
end

return

function [jx,jy,jz,e]=MF_moments_hyper(hvec,h_dipol,t,J,gLande,Hcf,A,I)
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
Ham=Hcfh+Hzeeman-h_dipol(1)*Jxh-h_dipol(2)*Jyh-h_dipol(3)*Jzh + A*(Ixh*Jxh + Iyh*Jyh + Izh*Jzh);

%Diagonalize
[v,e]=eig(Ham);
e=real(diag(e));
e=e-min(e);
[e,n]=sort(e);
v=v(:,n);

%Calculate Matrixelements and Moments
if t==0
    % At zero temperature, use only lowest eigenvalue.
    jx=real(v(:,1)'*Jxh*v(:,1));
    jy=real(v(:,1)'*Jyh*v(:,1));
    jz=real(v(:,1)'*Jzh*v(:,1));
else
    % Boltzman factor (with t in Kelvin)
    % energien korrigieren, damit positiv, sonst NaN Fehler mit exp()
    e=e-min(e);
    z=exp(-e/(t/11.6))/sum(exp(-e/(t/11.6)));
    jx=real(diag(v'*Jxh*v))'*z;
    jy=real(diag(v'*Jyh*v))'*z;
    jz=real(diag(v'*Jzh*v))'*z;
end

return
