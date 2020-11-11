function [E_0,Egs,Egs2,acc_rate,mag,lattice,lat_mom,tempMark,T,mark]=rand_walk(sweeps,pt_intv,ion,params,inter,lattice,T,E_0,lat_mom,field_change,numWkrs,tempMark,prob)

p=1; % measurement counters
mark = 2; % parallel tempering counters
kB = 1/11.6; % Boltzmann constant [meV/K]
N = params.N_Er; % Number of lattice sites

if iscell(lat_mom)
    lat_mom = lat_mom{:};
end

% persistent phasex phasey phase0 % local variables
% if isempty(phasex)
%     phase0=ones(bestE);
%     qx=2*pi/L;
%     qy=qx;
%     pos=zeros(3,bestE);
%     for i=1:bestE
%         pos(:,i)=lattice{i}.position;c
%     end
%     phasex=exp(1i*qx*pos(1,:))'*exp(1i*qx*pos(1,:));
%     phasey=exp(1i*qy*pos(2,:))'*exp(1i*qy*pos(2,:));
% end

% arrays to write the magnetizations
mag = zeros([3,4,sweeps]);
% magsq = zeros([3,3,4,N_meas]);

% % form factors
% sq0 = zeros([3,3,4,4,N_meas]);
% sqx = zeros([3,3,4,4,N_meas]);
% sqy = zeros([3,3,4,4,N_meas]);

% energies per spin and energies per spin squared
Egs = zeros([1,sweeps]);
Egs2 = zeros([1,sweeps]);
acc_rate = zeros([1,sweeps]);

for j=1:N
    lattice{j}.mom = lat_mom(1:3,j)';
    lattice{j}.beta = lat_mom(4,j);
    lattice{j}.alpha = lat_mom(5,j);
    lattice{j}.energy = lat_mom(6,j);
end

accept = 0; % Initiate acceptance count
for iterations=1:(sweeps*N)
    % Rotate one moment at a randomly picked site
    chosen = 0;
    while(chosen == 0)
        new_ion=randi(size(lattice,2),1,1);
        [alpha_new, beta_new] = random_angles;
        [new_mom, new_E_Zee,field_change] = cl_eff_mod(alpha_new, beta_new, params.field, lattice{new_ion},field_change);
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
    
    % take a measurements of <E>, <E^2>, and form factors every meas_intv steps
    if mod(iterations,N) == 0
        Egs(p) = (E_0/N);  %Energy per site
        Egs2(p) = (E_0/N)^2; %Energy squared per site
        acc_rate(p) = accept/N; % Calculate acceptance rate of the metropolis steps between measurements
        accept = 0;
        
        % magnetizations and magnetizations squared per type of site
        mag(:,1,p) = sum(lat_mom(1:3,1:4:end)/N,2);
        mag(:,2,p) = sum(lat_mom(1:3,2:4:end)/N,2);
        mag(:,3,p) = sum(lat_mom(1:3,3:4:end)/N,2);
        mag(:,4,p) = sum(lat_mom(1:3,4:4:end)/N,2);
        
%         % Structural factor
%         for m=1:4
%             for mp=1:4
%                 sq0(1:3,1:3,m,mp,p) = lat_mom(1:3,m:4:end).*phase0(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                 sqx(1:3,1:3,m,mp,p) = lat_mom(1:3,m:4:end).*phasex(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%                 sqy(1:3,1:3,m,mp,p) = lat_mom(1:3,m:4:end).*phasey(m:4:end,mp:4:end).*lat_mom(1:3,mp:4:end)';
%             end
%             magsq(1:3,1:3,m,p) = sum((lat_mom(1:3,m:4:end)./N).*(lat_mom(1:3,1:4:end)./bestE),2);
%         end

        p=p+1;
    end
%% Parallel tempering
    if pt_intv ~= inf % Determing wether or not to implement parallel tempering 
        token = tempMark(1);
        % Every pt_intv steps, perform a parallel tempering attempt (even and odd numbered workers' turns
        % are mismatched by half interval) To use this procedure, it is necessary to ensure all cores are working
        % at the same time throughtout (temperature points has to match the number of cores).
        % the energy comparison is done with global energy, in accordance with function metropolis(dE,T)
        if rem(iterations,pt_intv*N)==0
            if mod(mark,2) == 0
                fprintf('Parallel tempering attempt No.%u.\n',iterations/(pt_intv*N));
%                 fprintf('Parallel tempering attempt No.%u, current iteration: %u/%u.\n',mark,iterations,(pt_intv*N));
                % Even numbered workers compare and exchange with the nearst neighbour on the left
                if mod(labindex,2) == 0
                    left = mod(labindex-2,numWkrs)+1; % worker on the left
                    %                 fprintf('Exchanging with worker %u.\n',left); % For debugging
                    E_l = labSendReceive(left,left,E_0,1);
                    T_l = labSendReceive(left,left,T,2);
                    newToken = labSendReceive(left,left,token,3);
                    chance = exp( 1/kB*(E_0-E_l)*(1/T-1/T_l) );
                    fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
                    if prob(mark) <= chance
                        fprintf(1,'Swaping with the left neighbour.\n');
                        T = T_l;
                        tempMark(mark) = newToken; %Mark down the swap on MC timeline
                        token = newToken;
                    else
                        tempMark(mark) = token;
                        fprintf(2,'swap rejected.\n');
                    end
                    %             % Odd number workers compare and exchange with the nearst neighbour on the right
                else
                    right = mod(labindex,numWkrs)+1; % worker on the right
                    %                 fprintf('Exchanging with worker %u.\n',right); % For debugging
                    E_r = labSendReceive(right,right,E_0,1);
                    T_r = labSendReceive(right,right,T,2);
                    newToken = labSendReceive(right,right,token,3);
                    chance = exp( 1/kB*(E_0-E_r)*(1/T-1/T_r) );
                    fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
                    if prob(mark) <= chance
                        fprintf(1,'Swaping with the right neighbour.\n');
                        T = T_r;
                        token = newToken;
                    else
                        fprintf(2,'swap rejected.\n');
                    end
                end
                tempMark(mark) = token; %Mark down the swap on MC timeline
                mark = mark + 1;
            else
                % the same process but toward the opposite exchange direction
                fprintf('Parallel tempering attempt No.%u.\n',iterations/(pt_intv*N));
%                 fprintf('Parallel tempering attempt No.%u, current iteration: %u/%u.\n',mark,iterations,(pt_intv*N));
                % Even number workers compare and exchange with the nearst neighbour on the right
                if mod(labindex,2) == 0
                    right = mod(labindex,numWkrs)+1; % worker on the right
                    %                 fprintf('Exchanging with worker %u.\n',right); % For debugging
                    E_r = labSendReceive(right,right,E_0,1);
                    T_r = labSendReceive(right,right,T,2);
                    newToken = labSendReceive(right,right,token,3);
                    chance = exp( 1/kB*(E_0-E_r)*(1/T-1/T_r) );
                    fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
                    if prob(mark) <= chance
                        fprintf(1,'Swaping with the right neighbour.\n');
                        T = T_r;
                        token = newToken;
                    else
                        fprintf(2,'swap rejected.\n');
                    end
                    tempMark(mark) = token; %Mark down the swap on MC timeline
                    % Odd number workers compare and exchange with the nearst neighbour on the left
                else
                    left = mod(labindex-2,numWkrs)+1; % worker on the left
                    %                 fprintf('Exchanging with worker %u.\n',left); % For debugging
                    E_l = labSendReceive(left,left,E_0,1);
                    T_l = labSendReceive(left,left,T,2);
                    newToken = labSendReceive(left,left,token,3);
                    chance = exp( 1/kB*(E_0-E_l)*(1/T-1/T_l) );
                    fprintf('Criterion probability: %1$.3f, Transition probability: %2$.3f\n', prob(mark), chance); %checkpoint
                    if prob(mark) <= chance
                        fprintf(1,'Swaping with the left neighbour.\n');
                        T = T_l;
                        token = newToken;
                    else
                        fprintf(2,'swap rejected.\n');
                    end
                end
                tempMark(mark) = token; %Mark down the swap on MC timeline
                mark = mark + 1;
            end
        end
    end
end

% compress sampling points into single datapoints
Egs = mean(Egs); % Energy mean over all the statistical points
Egs2 = mean(Egs2); % Energy mean square over all the statistical points

clearvars -except E_0 Egs Egs2 acc_rate mag lattice lat_mom tempMark T mark p
return