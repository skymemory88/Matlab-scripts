% function [energies,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,E_0,lat_mom]=LiIonsF4_MC(ion,L,Niter,inter,lattice,field,T,E_0,lat_mom,field_change)
function [energies,Egs,Egs2,mag,lattice,lat_mom,tempMark,T]=LiIonsF4_MC(ion,params,inter,lattice,T,E_0,lat_mom,field_change,numWkrs,half_interval,tempMark,prob)

kB = 1/11.6; % Boltzmann constant
N_meas=floor(params.Nitermeas/params.N_Er);

%N_meas=1;
if iscell(lat_mom)
    lat_mom = lat_mom{:};
end

% persistent phasex phasey phase0 % local variables
% if isempty(phasex)
%     phase0=ones(params.N_Er);
%     qx=2*pi/L;
%     qy=qx;
%     pos=zeros(3,params.N_Er);
%     for i=1:params.N_Er
%         pos(:,i)=lattice{i}.position;
%     end
%     phasex=exp(1i*qx*pos(1,:))'*exp(1i*qx*pos(1,:));
%     phasey=exp(1i*qy*pos(2,:))'*exp(1i*qy*pos(2,:));
% end

energies=zeros(1,params.Nitermeas);

% arrays to write the magnetizations
mag=zeros([3,4,N_meas]);
% magsq=zeros([3,3,4,N_meas]);
% % form factors
% sq0=zeros([3,3,4,4,N_meas]);
% sqx=zeros([3,3,4,4,N_meas]);
% sqy=zeros([3,3,4,4,N_meas]);

% energies per spin and energies per spin squared
Egs=zeros([1,N_meas]);
Egs2=zeros([1,N_meas]);

% counters
p=1;
mark = 1;

for j=1:params.N_Er
    lattice{j}.mom=lat_mom(1:3,j)';
    lattice{j}.beta=lat_mom(4,j);
    lattice{j}.alpha=lat_mom(5,j);
    lattice{j}.energy=lat_mom(6,j);
end

for iterations=1:params.Nitermeas
    
    % Rotate one moment
    chosen = 0;
    while(chosen == 0)
        new_ion=randi(size(lattice,2),1,1);
        [alpha_new, beta_new]=random_angles;
        [new_mom, new_E_Zee,field_change]=cl_eff_mod(alpha_new, beta_new, params.field, lattice{new_ion},field_change);
        chosen = 1;
    end
    
    % Computes the energy variation
    dE = energy_update(ion,params.field,lattice,params.L,inter,new_ion,new_mom,new_E_Zee,lattice{new_ion}.energy);
    acc = metropolis(dE,T); % Metropolis-Hastings
    %fprintf('dE = %f meV, acc = %d \n',dE,acc)
    %acc=glauber(dE,T); % Glauber rate
    if acc == 1
        lattice{new_ion}.mom=new_mom;
        lattice{new_ion}.alpha=alpha_new;
        lattice{new_ion}.beta=beta_new;
        lattice{new_ion}.energy=new_E_Zee;
        lat_mom(1:3,new_ion)=new_mom;
        lat_mom(4,new_ion)=alpha_new;
        lat_mom(5,new_ion)=beta_new;
        lat_mom(6,new_ion)=new_E_Zee;
    end
    
    E_0 = E_0+acc*dE;
    energies(iterations) = E_0/params.N_Er;  %Energy per site
    
    % take a measure of <E>, <E^2>, and form factors every params.N_Er steps
    if rem(iterations,params.N_Er)==0
        Egs(p)=energies(iterations); %Energy per site
        Egs2(p)=(energies(iterations))^2; %Energy squared per site        
%         if p > 1
%             varE = mean(Egs2(1:p-1))-mean(Egs(1:p-1))^2;
%             fprintf('Variance of E: %.3f.\n',varE);
%             Ursell = (Egs(p-1)*Egs(p)) - mean(Egs(1:p-1))*mean(Egs(1:p));
%             fprintf('Ursell function: %.3f.\n',Ursell);
%             autoC = abs(Ursell/varE); % comput autocorrelation function
% %             fprintf('Autocorrelation function: %.3f.\n',autoC);
%         else
%             autoC = 1/exp(1);
%         end        
        %% Parallel tempering
        % Every 50*params.N_Er steps, perform a parallel tempering attempt (even and odd numbered workers' turns
        % are mismatched by half interval) To use this procedure, it is necessary to ensure all cores are working
        % at the same time throughtout (temperature points has to match the number of cores).
        if rem(iterations,2*half_interval)==0 && iterations > 0.5*params.Nitermeas % Only start parallel tempering after half the interations are done
            % Exchange to the left for even numbered worker
            if mod(labindex,2) == 0
                left = mod(labindex-2,numWkrs)+1; % worker on the left
                E_l = labSendReceive(left,left,Egs(p),1);
                T_l = labSendReceive(left,left,T,2);
                if prob(mark) <= exp( (Egs(p)-E_l)/kB/(T-T_l) )
                    fprintf('Swaping with %d.\n',left);
                    T = T_l;
                    tempMark(mark) = left; %Mark down the swap on MC timeline
                end
            else
                right = mod(labindex,numWkrs)+1; % worker on the right
                %                 [E_h,T_h] = labSendReceive(right,right,[Egs(p),T]);
                E_h = labSendReceive(right,right,Egs(p),1);
                T_h = labSendReceive(right,right,T,2);
                chance = exp( (Egs(p)-E_h)/kB/(T-T_h) );
                fprintf('Swaping probability: %1$.3f, Transition probability: %2$.3f', prob(mark), chance); %checkpoint
                if prob(mark) <= exp( (Egs(p)-E_h)/kB/(T-T_h) )
                    fprintf(', Swaping with %d.\n',right);
                    T = T_h;
                    tempMark(mark) = right; %Mark down the swap on MC timeline
                end
            end
            mark = mark + 1;
        elseif rem(iterations,2*half_interval+params.N_Er)==0
            % Exchange to the left for odd numbered worker
            if mod(labindex,2) ~= 0
                left = mod(labindex-2,numWkrs)+1; % worker on the left
                %                 [E_l,T_l] = labSendReceive(left,left,[Egs(p),T]);
                E_l = labSendReceive(left,left,Egs(p),1);
                T_l = labSendReceive(left,left,T,2);
                if prob(mark) <= exp( (Egs(p)-E_l)/kB/(T-T_l) )
                    fprintf('Swaping with %d.\n',left);
                    T = T_l;
                    tempMark(mark) = left; %Mark down the swap on MC timeline
                end
            else
                right = mod(labindex,numWkrs)+1; % worker on the right
                E_h = labSendReceive(right,right,Egs(p),1);
                T_h = labSendReceive(right,right,T,2);
                chance = exp( (Egs(p)-E_h)/kB/(T-T_h) );
                 fprintf('Swaping probability: %1$.3f, Transition probability: %2$.3f', prob(mark), chance); %checkpoint
                if prob(mark) <= exp( (Egs(p)-E_h)/kB/(T-T_h) )
                    fprintf(', Swaping with %d.\n',right);
                    T = T_h;
                    tempMark(mark) = right; %Mark down the swap on MC timeline
                end
            end
            mark = mark + 1;
        end
%         if autoC <= 1/exp(1)
            % magnetizations and magnetizations squared per type of site
            mag(:,1,p)=sum(lat_mom(1:3,1:4:end)/params.N_Er,2);
            mag(:,2,p)=sum(lat_mom(1:3,2:4:end)/params.N_Er,2);
            mag(:,3,p)=sum(lat_mom(1:3,3:4:end)/params.N_Er,2);
            mag(:,4,p)=sum(lat_mom(1:3,4:4:end)/params.N_Er,2);
            
%             % Structural factor
%             for m=1:4
%                 for mp=1:4
%                     sq0(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phase0(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                     sqx(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phasex(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                     sqy(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phasey(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                 end
%                 magsq(1:3,1:3,m,p)=sum((lat_mom(1:3,m:4:end)./params.N_Er).*(lat_mom(1:3,1:4:end)./params.N_Er),2);
%             end
      
            p=p+1;
%         end
    end
end
energies = energies(1:params.N_Er:end); % keep data only every N = Num of Spin steps
tempMark = tempMark(1:mark-1);
% clearvars -except energies Egs Egs2 mag magsq sq0 sqx sqy lattice E_0 lat_mom tempMark
clearvars -except energies Egs Egs2 mag lattice lat_mom tempMark T
return