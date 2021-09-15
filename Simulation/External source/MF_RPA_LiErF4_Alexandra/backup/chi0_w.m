function [chi0,e,v]=chi0_w(ion,hvec,t,momente,momente_hyp,ishf,dip_range,withdemagn,alpha,omega,epsilon)

%Initialize chi
chi0 = zeros(3,3,size(momente,1),length(omega));

%Initialize dipole and exchange terms (D in J*D*J)
q = [0 0 0];
[d_dip,d_ex] = calc_dip_ex(ion,q,dip_range,withdemagn,alpha);

%% Haven't standardized this between here and remf_3
for ionn=1:size(ion.tau,1)

    %calculate meanfield
     h_dipole=[0 0 0];
     h_ex=[0 0 0];
     h_dipole_hyp=[0 0 0];
     h_ex_hyp=[0 0 0];
     for ionm=1:size(ion.tau,1)
         h_dipole=h_dipole+momente(ionm,:)*diag(ion.renorm)*d_dip(:,:,ionm,ionn)';
         h_ex=h_ex+momente(ionm,:)*diag(ion.renorm)*d_ex(:,:,ionm,ionn)';
         h_dipole_hyp=h_dipole_hyp + momente_hyp(ionm,:)*diag(ion.renorm)*d_dip(:,:,ionm,ionn)';
         h_ex_hyp=h_ex_hyp + momente_hyp(ionm,:)*diag(ion.renorm)*d_ex(:,:,ionm,ionn)';
     end
     
     %h_mf is the virtual crystal mean field .
     h_mf=(1-ion.iso_HF)*(h_dipole+h_ex) + ion.iso_HF*(h_dipole_hyp+h_ex_hyp);
     
     %% Generates expectation values. 
     [mJx,mJy,mJz,e,ne,v] = gen_MF_moments_matrix(ion,hvec,h_mf,t,ishf); %% Done here without HF for Er, needs to be modified to be more general
 
 for n=1:length(omega)
    w=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)))./ ...
      (e*ones(1,size(v,1))-ones(size(v,1),1)*e'-omega(n)-1i*epsilon);

    chi0(:,:,ionn,n)=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
          sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
          sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];
 end  
end
return

