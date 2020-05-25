function [lattice,params,svec,lat_mom]=lat(ion,params)

L=params.L; % size of the cubic box
prop=params.prop; % proportions of each ions


[svec,lvec,kvec,hvec]=ndgrid(1:4, 1:L, 1:L, 1:L);
hkls=[hvec(:) kvec(:) lvec(:) svec(:)];
% hvec, kvec, lvec : position for each unit cell, svec : type of the site

params.N_Er=0;
params.N_Ho=0;
params.N_Yb=0;

lattice=cell(1,size(hkls,1));

pop_lat=randomlat(prop,size(hkls,1));  % lattice filling according to the ion proportions
lat_mom=zeros(6,size(hkls,1));

for i=1:size(hkls,1)
    lattice{i}=ion{pop_lat(i)}; % randomly fills the lattice in accordance with the proportions constraints
    lattice{i}.position=[hkls(i,1);hkls(i,2);hkls(i,3)]+params.pos(:,hkls(i,4));
    lattice{i}.site=hkls(i,4);
    lattice{i}.staggfield=params.staggfield(hkls(i,4),:);
    
    switch lattice{i}.name
        case 'Er'
            params.N_Er=params.N_Er+1;
            % initial configuration
            if(params.init_Er==2) % bi-layered AFM order
                if(lattice{i}.site==1 || lattice{i}.site==4)
                    [lattice{i}.mom,lattice{i}.energy,params.field_changed]=cl_eff_mod(pi/2, 0,[0.00001,0,0],lattice{i},params.field_changed);
                    lattice{i}.alpha=pi/2;
                    lattice{i}.beta=0;
                else
                    [lattice{i}.mom,lattice{i}.energy,params.field_changed]=cl_eff_mod(pi/2,pi,[0.00001,0,0],lattice{i},params.field_changed);
                    lattice{i}.alpha=pi/2;
                    lattice{i}.beta=pi;
                end
            end
            if(params.init_Er==1) % random orientations
                [lattice{i}.alpha, lattice{i}.beta]=random_angles;
                [lattice{i}.mom,lattice{i}.energy,params]=cl_eff_mod(lattice{i}.alpha,lattice{i}.beta, [0,0,0], lattice{i},params);
            end
            
        case 'Ho'
            params.N_Ho=params.N_Ho+1;
            % initial configuration
            if(params.init_Ho==2) % spins parallel to the z axis
                [lattice{i}.mom,lattice{i}.energy,params]=cl_eff_mod(pi, pi,[0,0,0],lattice{i},params);
                lattice{i}.alpha=pi;
                lattice{i}.beta=pi;
            end
            if(params.init_Ho==1) % random orientations
                [lattice{i}.alpha, lattice{i}.beta]=random_angles;
                [lattice{i}.mom,lattice{i}.energy, params]=cl_eff_mod(lattice{i}.alpha,lattice{i}.beta, [0,0,0], lattice{i},params);
            end
        case 'Yb'
            params.N_Yb=params.N_Yb+1;
            if(params.init_Yb==1) % random orientations
                [lattice{i}.alpha, lattice{i}.beta]=random_angles;
                [lattice{i}.mom,lattice{i}.energy,params]=cl_eff_mod(lattice{i}.alpha,lattice{i}.beta, [0,0,0], lattice{i},params);
            end
            if(params.init_Yb==2) % bi-layered AFM order
                if(lattice{i}.site==1 || lattice{i}.site==4)
                    [lattice{i}.mom,lattice{i}.energy,params]=cl_eff_mod(pi/2, 0,[0.00001,0,0],lattice{i},params);
                    lattice{i}.alpha=pi/2;
                    lattice{i}.beta=0;
                else
                    [lattice{i}.mom,lattice{i}.energy,params]=cl_eff_mod(pi/2,pi,[0.00001,0,0],lattice{i},params);
                    lattice{i}.alpha=pi/2;
                    lattice{i}.beta=pi;
                end
            end
            
    end
    lat_mom(1:3,i)=lattice{i}.mom;
    lat_mom(4,i)=lattice{i}.alpha;
    lat_mom(5,i)=lattice{i}.beta;
    lat_mom(6,i)=lattice{i}.energy;
end


