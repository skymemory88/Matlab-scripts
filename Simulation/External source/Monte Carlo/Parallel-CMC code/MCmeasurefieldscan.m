function MCmeasurefieldscan(filein_id, jobid)

load(sprintf('resultsEQ_%d.mat',filein_id));
   
% for f=1:size(params.field,1) % Appears to be duplicated loop
%    [bestE(:,f),bestE2(:,f),C_fdt(:,f),angl(f,:),mmagsq(:,:,:,f),msq0(:,:,:,:,f),msqx(:,:,:,:,f),msqy(:,:,:,:,f),mmag(:,:,f),malt(f),lat_mom,params]=Fieldloop(ion,params,inter,lattice,EQlat_mom{1,f},params.field(f,:));
% end

[bestE,bestE2,C_fdt,angl,mmagsq,msq0,msqx,msqy,mmag,malt,lat_mom,params]=paraFieldloop(ion,params,inter,lattice,EQlat_mom);

% DC susceptibility (derivative of the energy with respect to field) --Yikai
X0 = diff(bestE)./diff(params.field(:,3)');

% AC susceptibility 
chi_fdt=zeros(3,3,length(params.field));
chi=zeros(3,3,length(params.field));
for a=1:3
    for b=1:3
        for f=1:length(params.field)
            if(params.temp~=0)
%                 chi_fdt(a,b,f)=(mean(mmagsq(f,a,b,:))-mean(mmag(f,a,:))*mean(mmag(f,b,:)))/params.temp;
%                 chi(f,a,b)=mean(mean(msq0(f,a,b,:,:)))/params.temp;
                chi_fdt(a,b,f)=(mean(mmagsq(a,b,:,f))-mean(mmag(a,:,f))*mean(mmag(b,:,f)))/params.temp;
                chi(f,a,b)=mean(mean(msq0(a,b,:,:,f)))/params.temp;
            end
        end
    end
end
   cd('G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Quantum Monte Carlo\Tc scaling');
   filename=sprintf('results_fieldscan_%d_%i.mat',params.jobid,jobid);
   save(filename,'bestE','bestE2','C_fdt','X0','chi_fdt','chi','mmagsq','angl','msq0','msqx','msqy','mmag','malt','params','lat_mom');
   cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\External source\Monte Carlo\Parallel-CMC code');
end