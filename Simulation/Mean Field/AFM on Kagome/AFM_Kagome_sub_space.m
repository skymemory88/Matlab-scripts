sigx = [0 1; 1 0];
sigy = [0 -1i; 1i 0];
sigz = [1 0; 0 -1];

Sx_1 = sigx;
Sx_2 = 1/2.*sigx;

Sy_1 = sigy;
Sy_2 = 1/2.*sigy;

Sz_1 = sigz;
Sz_2 = 1/2.*sigz;

limit = 100; %set iteration limit
D = 1; %setting dimensions
k = 0.0862875; %Boltzmann constant

minField = -0.25;
maxField = 0.25;
fieldStep = 5e-4;

minTemp = 0.01;
maxTemp = 0.2;
tempStep = 0.04;

field = int8((maxField-minField)/fieldStep); % discretization of the field
% temperature = 0.05;
temperature = int8((maxTemp-minTemp)/tempStep); % discretization of the temperature
vector = zeros(int8(field*temperature),3);
% map = zeros(field, temperature);

delta = 0.001;
iterator = 1;
newAvSx_1 = 0;
newAvSx_2 = 0;

newAvSy_1 = 0;
newAvSy_2 = 0;

newAvSz_1 = 0;
newAvSz_2 = 0;
avm = zeros(temperature, field);
ti = 1;
for temp = minTemp:tempStep:maxTemp
    avSx_1 = 0; %initial guess of average Sx(1)
    avSx_2 = -0.5; %initial guess of average Sx(1/2)
    
    avSy_1 = 0; %initial guess of average Sy(1)
    avSy_2 = -0.5; %initial guess of average Sy(1/2)
    
    avSz_1 = 1.0; %initial guess of average Sz(1)
    avSz_2 = -0.5; %initial guess of average Sz(1/2)
    
    hi = 1;
    for hL = minField:fieldStep:maxField
        beta = 1/(k*temp);
        J_1 = 200/beta; %coupling strength
        J_2 = 2.0 * J_1;       
        while 1
            Hx = J_1*avSx_2*kron(Sx_2,Sx_1) + J_2*avSx_2*kron(Sx_2,Sx_2);
            Hy = J_1*avSy_2*kron(Sy_2,Sy_1) + J_2*avSy_2*kron(Sy_2,Sy_2);
            Hz = J_1*avSz_2*kron(Sz_2,Sz_1) + J_2*avSz_2*kron(Sz_2,Sz_2);
            
            Hamlt = Hx + Hy + Hz - hL*kron(eye(2),Sz_1) - 2 * hL*kron(eye(2),Sz_2);
%             Hamlt = J_1*avSz_2*kron(kron(Sz_2,Sz_2),Sz_1) + J_2*avSz_2*kron(kron(Sz_1,Sz_2),Sz_2) - hL*kron(eye(4),Sz_1) - 2 * hL*kron(eye(4),Sz_2);
%             Hamlt = 2 * J_1*avSz_2*kron(Sz_2,Sz_1) + J_2*avSz_2*kron(Sz_2,Sz_2) - hL*kron(eye(2),Sz_1) - 2 * hL*kron(eye(2),Sz_2);
            [v, E] = eig(Hamlt);
            Z = trace(exp(-beta.*E));
            
            newAvSx_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sx_1)*v)/Z;
            newAvSx_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sx_2)*v)/Z;
            diffx_1 = newAvSx_1 - avSx_1;
            diffx_2 = newAvSx_2 - avSx_2;
            
            newAvSy_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sy_1)*v)/Z;
            newAvSy_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sy_2)*v)/Z;
            diffy_1 = newAvSy_1 - avSy_1;
            diffy_2 = newAvSy_2 - avSy_2;
            
            newAvSz_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sz_1)*v)/Z;
            newAvSz_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(2),Sz_2)*v)/Z;
            diffz_1 = newAvSz_1 - avSz_1;
            diffz_2 = newAvSz_2 - avSz_2;
            
            if abs(diffx_1) < delta && abs(diffx_2) < delta && abs(diffy_1) < delta && abs(diffy_2) < delta && abs(diffz_1) < delta && abs(diffz_2) < delta || iterator > limit
                %             if abs(diffz_1) < delta && abs(diffz_2) < delta || iterator > limit
                break
            else
                avSx_1 = newAvSx_1;
                avSx_2 = newAvSx_2;
                
                avSy_1 = newAvSy_1;
                avSy_2 = newAvSy_2;
                
                avSz_1 = newAvSz_1;
                avSz_2 = newAvSz_2;
            end
        end
        avm(ti,hi) = newAvSz_1 + 2 * newAvSz_2;
        %         vector(iterator,:) = [temp, hL, newAvSz_1 + 2 * newAvSz_2];
        %         map(floor(iterator/(temperature+1))+1,rem(iterator,temperature+1)+1) = avm;
        iterator = iterator + 1;
        %         if mod(iterator, 10) == 0
        %            sprintf('Current iteration: %d, and Current temperature: %d', iterator, temp);
        %         end
        hi = hi + 1;
    end
    ti = ti + 1;
end
% avm = [nonzeros(avm(1,:)) nonzeros(avm(2,:))];
figure
hold on
axis([minField maxField -2.5 2.5]);
xlabel('External magnetic filed (H)');
ylabel('Magnetization per site');
for i = 1:temperature
    plot(minField:fieldStep:maxField, avm(i,:),'o-');
end
% plot(vector(:,2),vector(:,3));
% scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
% view([0 90]);
