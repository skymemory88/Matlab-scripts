function srpa = gen_scattering_crosssec_chi0(ion,chi0,qvec,omega,temp,dip_range,withdemagn,alpha)

for nq=1:size(qvec,1)
q=qvec(nq,:);
clear Skw
chi = chi_qw(ion,q,chi0,dip_range,withdemagn,alpha);
%chim=chi_qw(q,chi0m,Jex);

for nw=1:length(omega)
    if temp==0
        Skw(nw)=1/pi* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    else
        Skw(nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
            sum(sum((eye(3)-(q'*q)/(q*q')).* ...
                imag(chi(:,:,nw))));%-chim(:,:,nw,nq))));
    end
end

srpa(nq)=spec1d(omega,Skw,sqrt(abs(Skw)));
%srpa(nq)=cut(srpa(nq),[-0.1 0.15]);
end
