function [momente,momente_hyp,e,v,evolution]=remf(ion,hvec,t,momente,dip_range,withdemagn,alpha)

% set convergence criteria. Convergence is accepted when either:
% (Difference < ConvMin) or (Difference < ConvMax and Niter>NiterMin) or (Niter>NiterMax)
% This ensures very good convergence when no problem, and decent
% convergence when difficult.
ConvMax=1e-6;
ConvMin=1e-7;
NiterMin=1;
NiterMax=10000;

% Initialize moments and evolutions
momente_hyp=momente; %we will apply the virtual crystal approximation for the isotops of LiErF4
momente_mean=momente;
evolution.er(1,:,:)=momente;
evolution.er_hyp(1,:,:)=momente;

alt=[1  1  1
    -1 -1  1
    -1 -1  1
     1  1  1];
 
% Calculates dipole and exchange terms


% if(isempty(d_dip))
    q = [0 0 0];
    [d_dip,d_ex] = calc_dip_ex(ion,q,dip_range,withdemagn,alpha);
% end

for iterations=1:NiterMax
 % iterations 
 
 momente_old_hyp = momente_hyp;
 momente_old = momente_mean;
 momente_old_mean=(1-ion.iso_HF)*momente_old+ion.iso_HF*momente_old_hyp;

 for ionn=1:size(ion.tau,1)
     %calculate meanfield
     h_dipole=[0 0 0];
     h_ex=[0 0 0];
     h_dipole_hyp=[0 0 0];
     h_ex_hyp=[0 0 0];
     for ionm=1:size(ion.tau,1)
         h_dipole=h_dipole+momente_old(ionm,:)*diag(ion.renorm)*d_dip(:,:,ionm,ionn)';
         h_ex=h_ex+momente_old(ionm,:)*diag(ion.renorm)*d_ex(:,:,ionm,ionn)';
         h_dipole_hyp=h_dipole_hyp + momente_old_hyp(ionm,:)*diag(ion.renorm)*d_dip(:,:,ionm,ionn)';
         h_ex_hyp=h_ex_hyp + momente_old_hyp(ionm,:)*diag(ion.renorm)*d_ex(:,:,ionm,ionn)';
     end
     
     %h_mf is the virtual crystal mean field.
     h_mf=(1-ion.iso_HF)*(h_dipole+h_ex) + ion.iso_HF*(h_dipole_hyp+h_ex_hyp);

     %calculate moments of ions in a meanfield (ie diagonnalize the
     %hamiltonian)
     ishf = 1;
         if ion.iso_HF > 0
             [jx,jy,jz,e,~,v]=gen_MF_moments(ion,hvec,h_mf,t,ishf);
             momente_hyp(ionn,:)=update_moments([jx,jy,jz],evolution.er_hyp,ionn,iterations);
         end
         if ion.iso_HF ~= 1
             [jx,jy,jz,e,~,v]=gen_MF_moments(ion,hvec,h_mf,t,1-ishf);
             momente(ionn,:)=update_moments([jx,jy,jz],evolution.er,ionn,iterations);
         end
 end
 
  momente_mean=(1-ion.iso_HF)*momente+ion.iso_HF*momente_hyp;
  
  evolution.er(size(evolution.er,1)+1,:,:) = momente;
  evolution.er_hyp(size(evolution.er_hyp,1)+1,:,:) = momente_hyp;
  
 % Iterate untill convergence criterion reached
  Diff=sum(sum(abs(momente_old_mean-momente_mean)));
  if (Diff<ConvMin) || (Diff<ConvMax && iterations>NiterMin) 
    disp(num2str([t,hvec,iterations,mean(alt.*momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f)'))
    evolution.iterations=iterations;
    return
  end
 end

disp(num2str([t,hvec,iterations,mean(alt.*momente_mean)],'Temperature: %3.3f, Field: (%3.3f,%3.3f,%3.3f), Iterations: %3.3i, <(-1)^(i+j)J>=(%3.3f,%3.3f,%3.3f), fail to converge!'))
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
    J=Jnew;
%     Jold=squeeze(evolution(end-2*ndif+1-4,ionn,:))'; Why the additional "-4"? --Yikai (12.02.2020)
%     Jnow=squeeze(evolution(end-1*ndif+1-4,ionn,:))';
%     Jnew=squeeze(evolution(end-0*ndif+1-4,ionn,:))';
    Jold=squeeze(evolution(end-2*ndif+1,ionn,:))';
    Jnow=squeeze(evolution(end-1*ndif+1,ionn,:))';
    Jnew=squeeze(evolution(end-0*ndif,ionn,:))';
%    enNr=min(max(0,(Jnew-Jnow)./(Jnow-Jold)),0.998);
%    if isempty(find(enNr>=1,1))
    enNr=(Jnew-Jnow)./(Jnow-Jold);
    for nxyz=1:3 % for each of the three (dimensions) components of the Js --Yikai (10.02.2020)
       if enNr(nxyz)>0.01 && enNr(nxyz)<0.998
          feNNr=(Jnew-Jnow)./(enNr-1);
          J(nxyz)=J(nxyz)-feNNr(nxyz);
       end
    end

%    if isempty(find(enNr>=1,1))
%        feNNr=(Jnew-Jnow)./(enNr-1);
%        J=Jnow-feNNr;
%    else
%        J=Jnew;
%    end
else
    J=strategies.damping*squeeze(evolution(end,ionn,:))'+(1-strategies.damping)*Jnew;
end
return


