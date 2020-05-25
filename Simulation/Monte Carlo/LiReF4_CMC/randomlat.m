function pop_lat=randomlat(prop, N)
% function which returns a vector containing the number of the ion for 
% each site of the lattice (1 for Er, 2 for Ho, etc.). prop is the vector
% with the proportions and N the number of sites.

if(sum(prop)~=1)
    disp('The sum of the proportions must be 1 !')
    return
end

N_ions=length(prop); % number of ion types
pop_lat=zeros(N,1);

n=zeros(N_ions,1);
for type=1:N_ions
    n(type,1)=floor(N*prop(type));
end

diff=N-sum(n,1);
if(diff>0)
    for i=1:diff % diff < N_ions
        n(i)=n(i)+1;
    end
end
% disp(['diff corr =', num2str(N-sum(n))])
% disp(n)

N_empty=N;

for type=1:N_ions
    
    n_rem=n(type); % number of ions of the considered type that we still have to put in the lattice
    N_rem=N_empty; % number of sites remaining to visit for this type of ion
    proba=n_rem/N_empty; % probability to fill the considered site with this ion type.
    
    for site=1:N
        if(pop_lat(site)==0)% empty site, if non zero, we skip it, already filled.
            test=rand(1);
            if(test<proba) 
                pop_lat(site)=type;
                N_empty=N_empty-1;
                n_rem=n_rem-1;
            end
             N_rem=N_rem-1;
        end       
        if(N_rem>0)
            proba=n_rem/N_rem; % update the probability, in order to fill with the right number of ions
        end
    end
end