function E_0=energy0(ion,params,inter,lattice,EQlat_mom)
% computes the energy of the configuration at T = 0 K

mu_b=5.7883817555e-2; % Bohr magneton in meV/T
if iscell(EQlat_mom)
    EQlat_mom = EQlat_mom{:};
end

%arrays to write the energy of each ion and the sum of the energies of the
% dipolar interactions for each ion
E_ions=zeros(1,size(lattice,2));
sum_dipole=zeros(3,size(lattice,2));

for i=1:size(lattice,2)
    
    ion1=lattice{i};
    gLande1=ion{ion1.num}.gLande;
    staggered_field=ion1.staggfield;
    ion1.mom=gLande1*EQlat_mom(1:3,i)';
    ion1.energy=EQlat_mom(6,i);
    theta=atan2(ion1.mom(2),ion1.mom(1));
    F=ion1.F;
    
    for j=1:size(lattice,2)
        
        ion2=lattice{j};
        gLande2=ion{ion2.num}.gLande;
        ion2.mom=gLande2*EQlat_mom(1:3,j)';
        ion2.energy=EQlat_mom(6,j);

        r_ij=ion2.position-ion1.position;
        [idx]=latmod(params,r_ij); % find the index to the distance info stored in the array 'inter'
        D=inter(:,:,idx); % retrieve the corresponding dipolar interaction matrix
        %disp([i j])
        Dmom=D*ion2.mom';
        if(i==j)
            %disp(Dmom)
            Dmom=2*Dmom;
        end
        sum_dipole(:,i)=sum_dipole(:,i)+Dmom;
        %fprintf('i %d, j %d, D*mom(j) %f %f %f\n', i,j,Dmom(1),Dmom(2),Dmom(3))

    end
    %fprintf('ion %d sum_dipole %f %f %f \n', i,sum_dipole(1,i),sum_dipole(2,i),sum_dipole(3,i))
    Hstaggfield=-mu_b*(staggered_field(1)*ion1.mom(1)+staggered_field(2)*ion1.mom(2)+staggered_field(3)*ion1.mom(3));
    Hdipole=-ion1.mom*sum_dipole(:,i)/2; 
    H_corr=-F*cos(4*theta);
    E_ions(i)=Hdipole+ion1.energy+Hstaggfield+H_corr;
    %fprintf('ion %d Hdipole %f Hz+Hcf %f \n', i,Hdipole,ion1.energy)
end

E_0=sum(E_ions(:));

end