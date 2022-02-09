function [ion,history,E,V] = LiIonsF4(ion,temp,field,phi,theta,demagn,alpha)
%Compute the magnetization and the alternated magnetization for a range
%of temperatures and fields. (1-x) is the dilution. temp and h are the 
%temperature and field arrays. xHypIso is the proportion of nuclear moments
%isotops for Er. alpha is, in the case the sample is an ellipsiod (a,c,c),
%the ratio a/c (for demagnetization tensor).

global strategies rundipole Options;
t = tic;
    
% Matrix used to calculate order parameters from (alternating) moments
alt=[1  1  1
    -1  1  1
    -1 -1  1
     1 -1  1];  

ion_mom_init = ion.mom;
 
ion.Jmom=zeros([3,length(temp),size(field,2),size(ion.name,1)]);
ion.Jmom_norm=zeros([length(temp),size(field,2),size(ion.name,1)]);
ion.Jmom_hyp=zeros([3,length(temp),size(field,2),size(ion.name,1)]);
ion.Jmom_hyp_norm=zeros([length(temp),size(field,2),size(ion.name,1)]);

ion.altJmom=zeros([3,length(temp),size(field,2),size(ion.name,1)]);
ion.altJmom_hyp=zeros([3,length(temp),size(field,2),size(ion.name,1)]);

ion.Js=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);
ion.Js_hyp=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);

ion.altJs=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);
ion.altJs_hyp=zeros([4,3,length(temp),size(field,2),size(ion.name,1)]);

history=cell([length(temp),size(field,2)]);

idx = find(ion.prop);
if ion.hyp(idx) > 0
    H_dims = (2*ion.I(idx)+1)*(2*ion.J(idx)+1); % Claculate the dimension of the Hilbert space
elseif ion.hyp(idx) == 0
    H_dims = 2*ion.J(idx)+1; % Claculate the dimension of the Hilbert space
end
E = zeros(size(field,2),length(temp),H_dims);
V = zeros(size(field,2),length(temp),H_dims,H_dims);

switch Options.scantype % determine x-variable and data slicing direction
    case 'field'
        continu_var = vecnorm(field,2);
        discrt_var = temp;
    case 'temp'
        continu_var = temp;
        discrt_var = vecnorm(field,2);
end

for j = 1:length(discrt_var)
    % Reinitiate the angular moments at each temperature step (2020.01.02)
    if length(continu_var) > 1
        ion.mom = ion_mom_init; %Reinitialize the spin moments at each field
        ion.mom_hyp = ion.mom;
        rundipole = true;
    end
    for i = 1:length(continu_var)
        switch Options.scantype % determine x-variable and data slicing direction
            case 'field'
                Bidx = i;
                Tidx = j;
            case 'temp'
                Bidx = j;
                Tidx = i;
        end
        [ion,evolution,energy,v,~] = remf(ion,field(:,Bidx)',temp(Tidx),demagn,alpha);        
        for k = 1:size(ion.name,1)
            ion.altJmom(:,Tidx,Bidx,k) = mean(alt.*ion.mom(:,:,k));
            ion.altJmom_hyp(:,Tidx,Bidx,k) = mean(alt.*ion.mom_hyp(:,:,k));
            ion.Jmom(:,Tidx,Bidx,k) = mean(ion.mom(:,:,k));
            ion.Jmom_hyp(:,Tidx,Bidx,k) = mean(ion.mom_hyp(:,:,k)); 
%             ion.Jmom_norm(:,j,i,k) = norm(ion.mom(1,:,k)); % from originial code--Yikai (11.02.2020)
            ion.Jmom_norm(Tidx,Bidx,k) = norm(ion.Jmom(:,Tidx,Bidx,k));
            ion.Jmom_hyp_norm(Tidx,Bidx,k) = norm(ion.Jmom_hyp(:,Tidx,Bidx,k)); % addded on 11.02.2020 --Yikai
            ion.Js(:,:,Tidx,Bidx,k) = ion.mom(:,:,k);
            ion.Js_hyp(:,:,Tidx,Bidx,k) = ion.mom_hyp(:,:,k);
            ion.altJs(:,:,Tidx,Bidx,k) = alt.*ion.mom(:,:,k);
            ion.altJs_hyp(:,:,Tidx,Bidx,k) = alt.*ion.mom_hyp(:,:,k);        
            history{Tidx,Bidx} = evolution;
        end
        
        if strategies.powerlaw && Tidx>3 && Tidx<length(temp) % Guess starting values for next temperature using powerlaw
            beta=0.5; % assume exponent 1/2. Could generalize to use more points and flexible exponent
            inv_beta=1/beta;
            J1=zeros([4,3,size(ion.name,1)]);
            J2=zeros([4,3,size(ion.name,1)]);
            Tc=zeros([4,3,size(ion.name,1)]);
            Aa=zeros([4,3,size(ion.name,1)]);
            Jnew=zeros([4,3,size(ion.name,1)]);
            for k=1:size(ion.name,1)
                J1(:,:,k)=squeeze(history{Tidx,Bidx}.(ion.name{k})(end,:,:));
                J2(:,:,k)=squeeze(history{Tidx-1,Bidx}.(ion.name{k})(end,:,:));
                Tc(:,:,k)=((J1(:,:,k)./J2(:,:,k)).^inv_beta*temp(Tidx-1)-temp(Tidx))./((J1(:,:,k)./J2(:,:,k)).^inv_beta-1); % can generalize to use more points
                if isempty(find(((Tc(:,:,k)<=0) + (Tc(:,:,k) > temp(end))),1)) % only do this if Tc guess is reasonable
                    Aa(:,:,k)=J2(:,:,k).^inv_beta./(Tc(:,:,k)-temp(Tidx-1));
                    if isempty(find(~isfinite(Aa(:,:,k))+isnan(Aa(:,:,k)),1)) % problem if Tc happens to equal temp(j-1)
                        Jnew(:,:,k)=real((Aa(:,:,k).*(Tc(:,:,k)-temp(Tidx+1))).^beta).*(temp(Tidx+1)<Tc(:,:,k));
                        ion.mom(:,:,k)=Jnew(:,:,k);
                    end
                end
            end
        end
        E(Bidx,Tidx,:) = energy; % Eigen energy
        V(Bidx,Tidx,:,:) = v; % Eigen function
    end
% Save the data split by temperatures when they are multi-dimensional, otherwise save the data outside this function
    if length(discrt_var) >1 && length(continu_var) >1
        switch Options.scantype
            case 'field'
                eee = squeeze(E(:,j,:));
                vvv = squeeze(V(:,j,:,:));
                ttt = temp(Tidx);
                fff = field;
                tit=strcat('Hscan_Li',ion.name(ion.prop~=0), sprintf('F4_%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                    discrt_var(j), theta*180/pi, phi*180/pi,ion.hyp(idx)));
                save(fullfile(Options.filepath,char(tit)),'ttt','fff','eee','vvv','ion','-v7.3')
            case 'temp'
                eee = squeeze(E(j,:,:));
                vvv = squeeze(V(j,:,:,:));
                ttt = temp;
                fff = field(:,Bidx);
                tit=strcat('Tscan_Li',ion.name(ion.prop~=0), sprintf('F4_%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                    discrt_var(j), theta*180/pi, phi*180/pi,ion.hyp(idx)));
                save(fullfile(Options.filepath,char(tit)),'ttt','fff','eee','vvv','ion','-v7.3')
        end
    end
end

ion.Jmom = squeeze(ion.Jmom(:,:,:,:));
ion.Jmom_hyp = squeeze(ion.Jmom_hyp(:,:,:,:));

ion.Jmom_norm = squeeze(ion.Jmom_norm(:,:,:));
ion.Jmom_hyp_norm = squeeze(ion.Jmom_hyp_norm(:,:,:));

ion.altJmom = squeeze(ion.altJmom(:,:,:,:));
ion.altJmom_hyp = squeeze(ion.altJmom_hyp(:,:,:,:));

ion.Js = squeeze(ion.Js(:,:,:,:,:));
ion.Js_hyp = squeeze(ion.Js_hyp(:,:,:,:,:));

ion.altJs = squeeze(ion.altJs(:,:,:,:,:));
ion.altJs_hyp = squeeze(ion.altJs_hyp(:,:,:,:,:));

toc(t)