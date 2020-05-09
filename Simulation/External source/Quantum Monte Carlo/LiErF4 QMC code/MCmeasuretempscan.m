function MCmeasuretempscan(filein_id,jobid)

load(sprintf('resultsEQ_%d.mat',filein_id));

% for t=1:length(params.temp) % Appears to be duplicated loop
%     [bestE(t),bestE2(t),C_fdt(t),angl(t,:),mmagsq(t,:,:,:),msq0(t,:,:,:,:),msqx(t,:,:,:,:),msqy(t,:,:,:,:),mmag(t,:,:),malt(t),lat_mom,params]=Tloop(ion,params,inter,lattice,EQlat_mom{t,1},temp(t));
% end

[bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=Tloop(ion,params,inter,lattice,EQlat_mom);

% specific heat (derivative of the energy with respect to temperature) --Yikai
C_v = diff(bestE)./diff(params.temp);

% susceptibility 
chi_fdt=zeros(3,3,length(params.temp));
chi=zeros(3,3,length(params.temp));
for a=1:3
    for b=1:3
        for t=1:length(params.temp)
            if(params.temp(t)~=0)
                chi_fdt(a,b,t)=(mean(mmagsq(t,a,b,:))-mean(mmag(t,:,a))*mean(mmag(t,:,b)))/params.temp(t);
                chi(a,b,t)=mean(mean(msq0(t,a,b,:,:)))/params.temp(t);
            end
        end
    end
end

filename=sprintf('results_tempscan_%d_%i.mat',params.jobid,jobid);
save(filename,'bestE','bestE2','C_fdt','C_v','chi_fdt','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom');


end