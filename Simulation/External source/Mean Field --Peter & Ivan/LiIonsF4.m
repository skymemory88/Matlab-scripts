function [ion,history,E,V]=LiIonsF4(ion,temp,field,phi,demagn,alpha)
%Compute the magnetization and the alternated magnetization for a range
%of temperatures and fields. (1-x) is the dilution. temp and h are the 
%temperature and field arrays. xHypIso is the proportion of nuclear moments
%isotops for Er. alpha is, in the case the sample is an ellipsiod (a,c,c),
%the ratio a/c (for demagnetization tensor).

global strategies;
global rundipole;
global Options;

t=tic;
    
% Matrix used to calculate order parameters from (alternating) moments
alt=[1  1  1
    -1  1  1
    -1 -1  1
     1 -1  1];  

ion_mom_init = ion.mom;
 
ion.Jmom=zeros([3,length(temp),size(field,2),size(ion.name,1)]);
ion.Jmom_norm=zeros([length(temp),size(field,2),size(ion.name,1)]);
ion.Jmom_hyp=zeros([3,length(temp),size(field,2),size(ion.name,1)]);

ion.altJmom=zeros([3,length(temp),size(field,2),size(ion.name,1)]);
ion.altJmom_hyp=zeros([3,length(temp),size(field,2),size(ion.name,1)]);

ion.Js=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);
ion.Js_hyp=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);

ion.altJs=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);
ion.altJs_hyp=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);

history=cell([length(temp),size(field,2)]);
 
%%
% if length(temp) > size(field,2)
% for i = 1:size(field,2)
%     
%     if length(temp) > 1 %Ions' initial moments ion.mom(:,:,1)=[1 0 0; 1
% 0 0; -1 0 0; -1 0 0];      %Er % ion.mom(:,:,1)=[0 1 0; 0 1 0; 0 -1 0; 0
% -1 0];      %Er ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho %
% ion.mom(:,:,2)=[3.473 -0.045 0.1;  3.473 -0.045 0.1;  3.473 -0.045 0.1;
% 3.473 -0.045 0.1];  %Ho ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];
% Yb ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
% ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd ion.mom(:,:,6)=[0
% 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y ion.mom_hyp=ion.mom;
%         rundipole = true;
%     end
% 
%     for j = 1:length(temp)
%         [ion,evolution,energy,v]=remf(ion,field(:,i)',temp(j),demagn,alpha);
%         E(i,j,:)=energy; V(i,j,:,:)=v;
% 
% 
%         for k=1:size(ion.name,1)
%             ion.altJmom(:,j,i,k)=mean(alt.*ion.mom(:,:,k));
%             ion.altJmom_hyp(:,j,i,k)=mean(alt.*ion.mom_hyp(:,:,k));
%             ion.Jmom(:,j,i,k)=mean(ion.mom(:,:,k));
%             ion.Jmom_norm(j,i,k)=norm(ion.mom(1,:,k));
%             ion.Jmom_hyp(:,j,i,k)=mean(ion.mom_hyp(:,:,k));
%             ion.Js(:,:,j,i,k)=ion.mom(:,:,k);
%             ion.altJs(:,:,j,i,k)=alt.*ion.mom(:,:,k);
%             ion.altJs_hyp(:,:,j,i,k)=alt.*ion.mom_hyp(:,:,k);
%             
%             history{j,i}=evolution;
%         end
%         
%         if strategies.powerlaw && j>3 && j<length(temp) % Guess starting
%         values for next temperature using powerlaw
%             beta=0.5; % assume exponent 1/2. Could generalize to use more
%             points and flexible exponent ibeta=1/beta;
%             
%             J1=zeros([4,3,size(ion.name,1)]);
%             J2=zeros([4,3,size(ion.name,1)]);
%             Tc=zeros([4,3,size(ion.name,1)]);
%             Aa=zeros([4,3,size(ion.name,1)]);
%             Jnew=zeros([4,3,size(ion.name,1)]); for k=1:size(ion.name,1)
%                 J1(:,:,k)=squeeze(history{j,i}.(ion.name{k})(end,:,:));
%                 J2(:,:,k)=squeeze(history{j-1,i}.(ion.name{k})(end,:,:));
%                 Tc(:,:,k)=((J1(:,:,k)./J2(:,:,k)).^ibeta*temp(j-1)-temp(j))./((J1(:,:,k)./J2(:,:,k)).^ibeta-1);
%                 can generalize to use more points if
%                 isempty(find(((Tc(:,:,k)<=0) + (Tc(:,:,k) >
%                 temp(end))),1)) % only do this if Tc guess is reasonable
%                     Aa(:,:,k)=J2(:,:,k).^ibeta./(Tc(:,:,k)-temp(j-1)); if
%                     isempty(find(~isfinite(Aa(:,:,k))+isnan(Aa(:,:,k)),1))
%                     problem if Tc happens to be equal temp(j-1)
%                         Jnew(:,:,k)=real((Aa(:,:,k).*(Tc(:,:,k)-temp(j+1))).^beta).*(temp(j+1)<Tc(:,:,k));
%                         ion.mom(:,:,k)=Jnew(:,:,k);
%                     end
%                 end
%             end
%         end
%     end
%         cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\saved_mf_results\3.42GHz\absorption')
%         save('fieldscan22.mat','temp','field','ion','history','E','V','-v7.3')
%         cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature')
% end
% 
% else % length(temp) < size(field,2)
%%
for j = 1:length(temp)
    % Reinitiate the angular moments at each temperature step (2020.01.02)
    if size(field,2) > 1
        ion.mom = ion_mom_init; %Reinitialize the spin moments at each field
        ion.mom_hyp = ion.mom;
        rundipole = true;
    end

    for i = 1:size(field,2)
        [ion,evolution,energy,v,h_mf1]=remf(ion,field(:,i)',temp(j),demagn,alpha);
        E(i,j,:)=energy;
        V(i,j,:,:)=v;
        h_mf2(i,:)=h_mf1(2,:);

        for k=1:size(ion.name,1)
            ion.altJmom(:,j,i,k)=mean(alt.*ion.mom(:,:,k));
            ion.altJmom_hyp(:,j,i,k)=mean(alt.*ion.mom_hyp(:,:,k));
            ion.Jmom(:,j,i,k)=mean(ion.mom(:,:,k));
            ion.Jmom_hyp(:,j,i,k)=mean(ion.mom_hyp(:,:,k)); 
%             ion.Jmom_norm(:,j,i,k)=norm(ion.mom(1,:,k)); % from originial code--Yikai (11.02.2020)
            ion.Jmom_norm(j,i,k)=norm(ion.Jmom(:,j,i,k));
            ion.Jmom_hyp_norm(j,i,k)=norm(ion.Jmom_hyp(:,j,i,k)); % addded on 11.02.2020 --Yikai
            ion.Js(:,:,j,i,k)=ion.mom(:,:,k);
            ion.Js_hyp(:,:,j,i,k)=ion.mom_hyp(:,:,k);
            ion.altJs(:,:,j,i,k)=alt.*ion.mom(:,:,k);
            ion.altJs_hyp(:,:,j,i,k)=alt.*ion.mom_hyp(:,:,k);
            
            history{j,i}=evolution;
        end
        
        if strategies.powerlaw && j>3 && j<length(temp) % Guess starting values for next temperature using powerlaw
            beta=0.5; % assume exponent 1/2. Could generalize to use more points and flexible exponent
            ibeta=1/beta;
            
            J1=zeros([4,3,size(ion.name,1)]);
            J2=zeros([4,3,size(ion.name,1)]);
            Tc=zeros([4,3,size(ion.name,1)]);
            Aa=zeros([4,3,size(ion.name,1)]);
            Jnew=zeros([4,3,size(ion.name,1)]);
            for k=1:size(ion.name,1)
                J1(:,:,k)=squeeze(history{j,i}.(ion.name{k})(end,:,:));
                J2(:,:,k)=squeeze(history{j-1,i}.(ion.name{k})(end,:,:));
                Tc(:,:,k)=((J1(:,:,k)./J2(:,:,k)).^ibeta*temp(j-1)-temp(j))./((J1(:,:,k)./J2(:,:,k)).^ibeta-1); % can generalize to use more points
                if isempty(find(((Tc(:,:,k)<=0) + (Tc(:,:,k) > temp(end))),1)) % only do this if Tc guess is reasonable
                    Aa(:,:,k)=J2(:,:,k).^ibeta./(Tc(:,:,k)-temp(j-1));
                    if isempty(find(~isfinite(Aa(:,:,k))+isnan(Aa(:,:,k)),1)) % problem if Tc happens to be equal temp(j-1)
                        Jnew(:,:,k)=real((Aa(:,:,k).*(Tc(:,:,k)-temp(j+1))).^beta).*(temp(j+1)<Tc(:,:,k));
                        ion.mom(:,:,k)=Jnew(:,:,k);
                    end
                end
            end
        end
        eee(i,1,:)=energy;
        vvv(i,1,:,:)=v;
    end
    ttt = temp(j);
% Save the data split by temperatures when there are multi-dimensional, otherwise save data outside this function
    if Options.saving == true
        if size(field,2) >1 && length(temp) >1
            cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4\output')
            tit=strcat('Hscan_Li',[ion.name(ion.prop~=0)], 'F4_', sprintf('%1$3.3fK_%2$uDeg',temp(j),phi*pi/180),'.mat');
            save(char(tit),'ttt','field','eee','vvv','h_mf2')
%             save(char(tit),'ttt','field','eee','vvv','h_mf2','-v7.3')
            cd('G:\My Drive\File sharing\Programming scripts\Matlab\Simulation\Mean Field\LiReF4')
        end
    end
end

% end

for k=1:size(ion.name,1)    
    ion.Jmom(:,:,:,k)=squeeze(ion.Jmom(:,:,:,k));
    ion.Jmom_hyp(:,:,:,k)=squeeze(ion.Jmom_hyp(:,:,:,k));
    
    ion.Jmom_norm(:,:,k)=squeeze(ion.Jmom_norm(:,:,k));
    ion.Jmom_hyp_norm(:,:,k)=squeeze(ion.Jmom_hyp_norm(:,:,k));
    
    ion.altJmom(:,:,:,k)=squeeze(ion.altJmom(:,:,:,k));
    ion.altJmom_hyp(:,:,:,k)=squeeze(ion.altJmom_hyp(:,:,:,k));
    
    ion.Js(:,:,:,:,k)=squeeze(ion.Js(:,:,:,:,k));
    ion.Js_hyp(:,:,:,:,k)=squeeze(ion.Js_hyp(:,:,:,:,k));
    
    ion.altJs(:,:,:,:,k)=squeeze(ion.altJs(:,:,:,:,k));
    ion.altJs_hyp(:,:,:,:,k)=squeeze(ion.altJs_hyp(:,:,:,:,k));
end


toc(t)