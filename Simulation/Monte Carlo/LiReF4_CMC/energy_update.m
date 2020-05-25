function dE=energy_update(ion,field,lattice,system_size,inter,new_ion,new_mom, new_E_Zee, old_E_Zee)
% update to the energy
% ones has to make sure that the effective moment has been computed and
% saved in lattice{i} for each site and that new_mom is the effective new
% moment, old_E_Zee and new_E_Zee are the energies computed by cl_eff_mod