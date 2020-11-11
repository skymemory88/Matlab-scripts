% function [energies,Egs,Egs2,mag,magsq,sq0,sqx,sqy,lattice,E_0,lat_mom]=LiIonsF4_MC(ion,L,Niter,inter,lattice,field,T,E_0,lat_mom,field_change)
function [Egs,Egs2,acc_rate,mag,lattice,lat_mom,tempMark,T,mark,p]=LiIonsF4_MC(ion,params,inter,lattice,T,E_0,lat_mom,field_change,numWkrs,tempMark,prob)

% kB = 1/11.6; % Boltzmann constant
N = params.N_Er;
N_meas = params.Nitermeas;
iterMax = N_meas * params.meas_intv * N;
energies=zeros(1,N_meas);% array for the energy per spin
% token = tempMark(1);
%N_meas=1; % for debugging

if iscell(lat_mom)
    lat_mom = lat_mom{:};
end

% persistent phasex phasey phase0 % local variables
% if isempty(phasex)
%     phase0=ones(N);
%     qx=2*pi/L;
%     qy=qx;
%     pos=zeros(3,N);
%     for i=1:N
%         pos(:,i)=lattice{i}.position;c
%     end
%     phasex=exp(1i*qx*pos(1,:))'*exp(1i*qx*pos(1,:));
%     phasey=exp(1i*qy*pos(2,:))'*exp(1i*qy*pos(2,:));
% end

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
acc_rate=zeros([1,N_meas]);

% counters
p=1;
mark = 2;

for j=1:N
    lattice{j}.mom=lat_mom(1:3,j)';
    lattice{j}.beta=lat_mom(4,j);
    lattice{j}.alpha=lat_mom(5,j);
    lattice{j}.energy=lat_mom(6,j);
end

accept = 0; % Initiate acceptance count
for iterations=1:iterMax
    % Rotate one moment at a randomly picked site
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
        accept = accept + 1;
    end
    
    E_0 = E_0+acc*dE;
    energies(mod(iterations-1,params.meas_intv*N)+1) = E_0/N;  %Energy per site
    
    % take a measure of <E>, <E^2>, and form factors every meas_intv steps
    if rem(iterations, params.meas_intv*N) == 0
        Egs(p)=energies(mod(iterations-1,params.meas_intv*N)+1); %Energy per site
        Egs2(p)=(energies(mod(iterations-1,params.meas_intv*N)+1))^2; %Energy squared per site 
        acc_rate(p) = accept/(params.meas_intv*N); % Calculate acceptance rate of the metropolis steps between measurements
        accept = 0; % reset acceptance count to zero
%% Parallel tempering
%         % Every params.pt_intv steps, perform a parallel tempering attempt (even and odd numbered workers' turns
%         % are mismatched by half interval) To use this procedure, it is necessary to ensure all cores are working
%         % at the same time throughtout (temperature points has to match the number of cores).
%         if rem(iterations,params.pt_intv*params.meas_intv*N)==0
%             fprintf('Parallel tempering attempt No.%u.\n',iterations/(params.pt_intv*params.meas_intv*N));
% %             % Even numbered workers compare and exchange with the nearst neighbour on the left
%             if mod(labindex,2) == 0
%                 left = mod(labindex-2,numWkrs)+1; % worker on the left
% %                 fprintf('Exchanging with worker %u.\n',left); % For debugging
%                 E_l = labSendReceive(left,left,Egs(p),1);
%                 T_l = labSendReceive(left,left,T,2);
%                 newToken = labSendReceive(left,left,token,3);
%                 chance = exp( (Egs(p)-E_l)*(1/kB/T-1/kB/T_l) );
%                 fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
%                 if prob(mark) <= exp( (Egs(p)-E_l)*(1/kB/T-1/kB/T_l) )  
%                     fprintf(', Swaping with %d.\n',left);
%                     T = T_l;
%                     tempMark(mark) = newToken; %Mark down the swap on MC timeline
%                     token = newToken;
%                 else
%                     tempMark(mark) = token;
%                     fprintf(2,', swap rejected.\n');
%                 end
% %             % Odd number workers compare and exchange with the nearst neighbour on the right
%             else
%                 right = mod(labindex,numWkrs)+1; % worker on the right
% %                 fprintf('Exchanging with worker %u.\n',right); % For debugging
%                 E_r = labSendReceive(right,right,Egs(p),1);
%                 T_r = labSendReceive(right,right,T,2);
%                 newToken = labSendReceive(right,right,token,3);
%                 chance = exp( (Egs(p)-E_r)*(1/kB/T-1/kB/T_r) );
%                 fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
%                 if prob(mark) <= exp( (Egs(p)-E_r)*(1/kB/T-1/kB/T_r) )
%                     fprintf(', Swaping with %d.\n',right);
%                     T = T_r;
%                     tempMark(mark) = newToken; %Mark down the swap on MC timeline
%                     token = newToken;
%                 else
%                     tempMark(mark) = token;
%                     fprintf(2,', swap rejected.\n');
%                 end
%             end
%             % Repeat the process but toward the opposite exchange direction
%             % Even number workers compare and exchange with the nearst neighbour on the right
%             if mod(labindex,2) == 0
%                 right = mod(labindex,numWkrs)+1; % worker on the right
% %                 fprintf('Exchanging with worker %u.\n',right); % For debugging
%                 E_r = labSendReceive(right,right,Egs(p),1);
%                 T_r = labSendReceive(right,right,T,2);
%                 newToken = labSendReceive(right,right,token,3);
%                 chance = exp( (Egs(p)-E_r)*(1/kB/T-1/kB/T_r) );
%                 fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
%                 if prob(mark) <= exp( (Egs(p)-E_r)*(1/kB/T-1/kB/T_r) )
%                     fprintf(', Swaping with %d.\n',right);
%                     T = T_r;
%                     tempMark(mark) = newToken; %Mark down the swap on MC timeline
%                     token = newToken;
%                 else
%                     tempMark(mark) = token;
%                     fprintf(2,', swap rejected.\n');
%                 end
%             % Odd number workers compare and exchange with the nearst neighbour on the left
%             else
%                 left = mod(labindex-2,numWkrs)+1; % worker on the left
% %                 fprintf('Exchanging with worker %u.\n',left); % For debugging
%                 E_l = labSendReceive(left,left,Egs(p),1);
%                 T_l = labSendReceive(left,left,T,2);
%                 newToken = labSendReceive(left,left,token,3);
%                 chance = exp( (Egs(p)-E_l)*(1/kB/T-1/kB/T_l) );
%                 fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
%                 if prob(mark) <= exp( (Egs(p)-E_l)*(1/kB/T-1/kB/T_l) )  
%                     fprintf(', Swaping with %d.\n',left);
%                     T = T_l;
%                     tempMark(mark) = newToken; %Mark down the swap on MC timeline
%                     token = newToken;
%                 else
%                     tempMark(mark) = token;
%                     fprintf(2,', swap rejected.\n');
%                 end
%             end
%             mark = mark + 1;
%         end
%%
        % magnetizations and magnetizations squared per type of site
        mag(:,1,p)=sum(lat_mom(1:3,1:4:end)/N,2);
        mag(:,2,p)=sum(lat_mom(1:3,2:4:end)/N,2);
        mag(:,3,p)=sum(lat_mom(1:3,3:4:end)/N,2);
        mag(:,4,p)=sum(lat_mom(1:3,4:4:end)/N,2);

%             % Structural factor
%             for m=1:4
%                 for mp=1:4
%                     sq0(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phase0(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                     sqx(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phasex(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                     sqy(1:3,1:3,m,mp,p)=lat_mom(1:3,m:4:end).*phasey(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                 end
%                 magsq(1:3,1:3,m,p)=sum((lat_mom(1:3,m:4:end)./N).*(lat_mom(1:3,1:4:end)./N),2);
%             end

        p=p+1;
    end
end
% clearvars -except energies Egs Egs2 mag magsq sq0 sqx sqy lattice E_0 lat_mom tempMark T mark p
clearvars -except Egs Egs2 acc_rate mag lattice lat_mom tempMark T mark p
return