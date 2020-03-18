%% calculate energy without double counting
function [E,S] = calculate_energy(S,h,J,spin,dim)

J1  = J(1);           % J1z/J1z
J2  = J(2)./J(1);     % J2z/J1z
J3  = J(3)./J(1);     % J2z/J1z

spin2 = spin*spin; 
S = [S; S(1,:)]; % enlarge due to periodic boundary conditions
S = [S S(:,1)];
E = 0;

for n=1:dim(1)
    for m=1:dim(2)
        E = E - h*spin*S(n,m);
        E = E + J1*spin2*S(n,m)*S(n,m+1)+J1*spin2*S(n,m)*S(n+1,m);
        E = E + J2*spin2*S(n,m)*S(n+1,m+1);
        E = E + J3*spin2*S(n+1,m)*S(n,m+1);
    end
end

E = E./(dim(1)*dim(2));
end
