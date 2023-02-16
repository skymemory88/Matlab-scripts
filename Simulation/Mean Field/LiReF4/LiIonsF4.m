function [ion,history,E,V] = LiIonsF4(ion,temp,field,phi,theta,alpha)
%Compute the magnetization and the alternated magnetization for a range
%of temperatures and fields. (1-x) is the dilution. temp and h are the 
%temperature and field arrays. xHypIso is the proportion of nuclear moments
%isotops for Er. alpha is, in the case the sample is an ellipsiod (a,c,c),
%the ratio a/c (for demagnetization tensor).

global strategies rundipole Options;
t = tic;
    
% Matrix used to calculate order parameters from (alternating) moments
alt = [1  1  1
      -1  1  1
      -1 -1  1
       1 -1  1];  

ion_mom_init = ion.mom;
ion.demag = Options.demag;
history = cell([length(temp),size(field,2)]);

idx = find(ion.prop);
if ion.hyp(idx) > 0
    H_dims = (2*ion.I(idx)+1)*(2*ion.J(idx)+1); % Calculate the dimension of the Hilbert space
elseif ion.hyp(idx) == 0
    H_dims = 2*ion.J(idx)+1; % Claculate the dimension of the Hilbert space
end

switch Options.scanMode % determine x-variable and data slicing direction
    case 'field'
        continu_var = vecnorm(field,2);
        discrt_var = temp;
    case 'temp'
        continu_var = temp;
        discrt_var = vecnorm(field,2);
end
E = zeros(length(continu_var),H_dims);
V = zeros(length(continu_var),H_dims,H_dims);

for ii = 1:length(discrt_var)
    % Reinitiate the angular moments at each temperature step (2020.01.02)
    if length(continu_var) > 1
        ion.mom = ion_mom_init; %Reinitialize the spin moments at each field
        ion.mom_hyp = ion.mom;
        rundipole = true;
    end
    ion.Jmom = zeros([3,length(continu_var),size(ion.name,1)]);
    ion.Jmom_norm = zeros([length(continu_var),size(ion.name,1)]);
    ion.Jmom_hyp = zeros([3,length(continu_var),size(ion.name,1)]);
    ion.Jmom_hyp_norm = zeros([length(continu_var),size(ion.name,1)]);

    ion.altJmom = zeros([3,length(continu_var),size(ion.name,1)]);
    ion.altJmom_hyp = zeros([3,length(continu_var),size(ion.name,1)]);

    ion.Js = zeros([4,3,length(continu_var),size(ion.name,1)]);
    ion.Js_hyp = zeros([4,3,length(continu_var),size(ion.name,1)]);

    ion.altJs = zeros([4,3,length(continu_var),size(ion.name,1)]);
    ion.altJs_hyp = zeros([4,3,length(continu_var),size(ion.name,1)]);
    for jj = 1:length(continu_var)
        switch Options.scanMode % determine x-variable and data slicing direction
            case 'field'
                [ion, evolution, energy, v, ~] = remf(ion, field(:,jj)', temp(ii), alpha);
            case 'temp'
                [ion, evolution, energy, v, ~] = remf(ion, field(:,ii)', temp(jj), alpha);
        end    
        for kk = 1:size(ion.name,1)
            ion.altJmom(:,jj,kk) = mean(alt.*ion.mom(:,:,kk));
            ion.altJmom_hyp(:,jj,kk) = mean(alt.*ion.mom_hyp(:,:,kk));
            ion.Jmom(:,jj,kk) = mean(ion.mom(:,:,kk));
            ion.Jmom_hyp(:,jj,kk) = mean(ion.mom_hyp(:,:,kk)); 
%             ion.Jmom_norm(:,j,i,k) = norm(ion.mom(1,:,k)); % from originial code--Yikai (11.02.2020)
            ion.Jmom_norm(jj,kk) = norm(ion.Jmom(:,jj,kk));
            ion.Jmom_hyp_norm(jj,kk) = norm(ion.Jmom_hyp(:,jj,kk)); % addded on 11.02.2020 --Yikai
            ion.Js(:,:,jj,kk) = ion.mom(:,:,kk);
            ion.Js_hyp(:,:,jj,kk) = ion.mom_hyp(:,:,kk);
            ion.altJs(:,:,jj,kk) = alt.*ion.mom(:,:,kk);
            ion.altJs_hyp(:,:,jj,kk) = alt.*ion.mom_hyp(:,:,kk);        
            history{ii,jj} = evolution;
        end
        if strcmp(Options.scanMode, 'field') % Guess starting values for next temperature using powerlaw
            if strategies.powerlaw && ii < length(discrt_var)
                beta = 0.5; % assume exponent 1/2. Could generalize to use more points and flexible exponent
                inv_beta = 1/beta;
                J1 = zeros([4,3,size(ion.name,1)]);
                J2 = zeros([4,3,size(ion.name,1)]);
                Tc = zeros([4,3,size(ion.name,1)]);
                Aa = zeros([4,3,size(ion.name,1)]);
                Jnew = zeros([4,3,size(ion.name,1)]);
                for kk = 1:size(ion.name,1)
                    J1(:,:,kk) = squeeze(history{ii,jj}.(ion.name{kk})(end,:,:));
                    J2(:,:,kk) = squeeze(history{ii-1,jj}.(ion.name{kk})(end,:,:));
                    Tc(:,:,kk) = ((J1(:,:,kk)./J2(:,:,kk)).^inv_beta*temp(ii-1)-temp(ii))./((J1(:,:,kk)./J2(:,:,kk)).^inv_beta-1); % can generalize to use more points
                    if isempty(find(((Tc(:,:,kk) <= 0) + (Tc(:,:,kk) > temp(end))),1)) % only do this if Tc guess is reasonable
                        Aa(:,:,kk) = J2(:,:,kk).^inv_beta./(Tc(:,:,kk)-temp(ii-1));
                        if isempty(find(~isfinite(Aa(:,:,kk))+isnan(Aa(:,:,kk)),1)) % problem if Tc happens to equal temp(j-1)
                            Jnew(:,:,kk) = real((Aa(:,:,kk).*(Tc(:,:,kk)-temp(ii+1))).^beta).*(temp(ii+1)<Tc(:,:,kk));
                            ion.mom(:,:,kk) = Jnew(:,:,kk);
                        end
                    end
                end
            end
        end
        E(jj,:) = energy; % Eigen energy
        V(jj,:,:) = v; % Eigen function
    end
    eee = squeeze(E(:,:));
    vvv = squeeze(V(:,:,:));
% Save the data split by temperatures when they are multi-dimensional, otherwise save the data outside this function
    if length(discrt_var) > 1 && length(continu_var) > 1
        switch Options.scanMode
            case 'field'
                ttt = temp(ii);
                fff = field;
                tit = strcat('Hscan_Li',ion.name(ion.prop ~= 0), sprintf('F4_%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                    discrt_var(ii), theta*180/pi, phi*180/pi,ion.hyp(idx)));
                save(fullfile(Options.filepath,char(tit)),'ttt','fff','eee','vvv','ion','-v7.3')
            case 'temp'
                ttt = temp;
                fff = field(:,ii);
                tit = strcat('Tscan_Li',ion.name(ion.prop ~= 0), sprintf('F4_%1$3.3fT_%2$.2fDg_%3$.1fDg_hp=%4$.2f.mat',...
                    discrt_var(ii), theta*180/pi, phi*180/pi,ion.hyp(idx)));
                save(fullfile(Options.filepath,char(tit)),'ttt','fff','eee','vvv','ion','-v7.3')
        end
    end
end

toc(t)