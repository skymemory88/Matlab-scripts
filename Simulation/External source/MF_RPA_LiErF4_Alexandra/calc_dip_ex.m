function [d_dip,d_ex] = calc_dip_ex(ion,q,dip_range,withdemagn,alpha)

% Dipole
eins=zeros(3,3,size(ion.tau,1),size(ion.tau,1)); eins(1,1,:,:)=1; eins(2,2,:,:)=1; eins(3,3,:,:)=1;
demagn=zeros(3,3,size(ion.tau,1),size(ion.tau,1));
demagn_t=ellipsoid_demagn(alpha);
demagn(1,1,:,:)=demagn_t(1,1);
demagn(2,2,:,:)=demagn_t(2,2);
demagn(3,3,:,:)=demagn_t(3,3);

%C = [3.48 3.48 0.94]; %% From Thesis C. Kramer eq. 4.4.1 (I don't think I need this, as without it the critical field is correct to what's reported

d_dip=ion.gLande^2*(0.05368*dipole_direct(ion,q,dip_range)+eins*ion.Lorenz/4-withdemagn*0.000745836*pi*demagn);
d_ex=exchange(ion,q);

return