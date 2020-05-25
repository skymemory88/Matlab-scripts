clear
format long g;

sigx = [0 1; 1 0];
sigy = [0 1i; -1i 0];
sigz = [1 0; 0 -1];

Sx = 1/2.*sigx;
Sy = 1/2.*sigy;
Sz = 1/2.*sigz;

Ix = 1/2.*sigx;
Iy = 1/2.*sigy;
Iz = 1/2.*sigz;

%N = 100; %setting site number
D = 1; %setting dimensions
k = 0.0862875;
J = 1.0; %coupling strength

minField = 0.01;
maxField = 0.6;
fieldStep = 0.004;

minTemp = 0.01;
maxTemp = 3;
tempStep = 0.02;

field = int8((maxField-minField)/fieldStep); % discretization of the field
temperature = int8((maxTemp-minTemp)/tempStep); % discretization of the temperature
vector = zeros(int8(field*temperature),3);
map = zeros(field, temperature);

delta = 0.001;
iterator = 1;
A = 0.0;
newAvSz = 0;

for hT = minField:fieldStep:maxField
    for temp = minTemp:tempStep:maxTemp
        avSz = 0.1; %initial guess of average Sz
        while 1
            Hamlt = -J*avSz*kron(eye(2),Sz) + hT*kron(eye(2),Sx) + A.*(kron(Ix,Sx)+kron(Iy,Sy)+kron(Iz,Sz));
            [v, E] = eig(Hamlt);
            beta = 1/(k*temp);
            Z = trace(exp(-beta.*E)); %Calculate the partition function weight
            newAvSz = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sz)*v)/Z; %Direct product by two because of the two neighbours
            if abs(newAvSz - avSz) < delta
                break
            else
                avSz = newAvSz;
            end
        end
        vector(iterator,:) = [temp, hT, avSz];
        map(floor(iterator/(temperature+1))+1,rem(iterator,temperature+1)+1) = avSz;
        iterator = iterator + 1;
        %if mod(iterater, 1000) == 0;
        %    iterater, temp, hT
        %end
    end
end

scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled')
view([0 90]);
