clear
format long g;

sigx = [0 1; 1 0];
sigy = [0 -1i; 1i 0];
sigz = [1 0; 0 -1];

sigx3 = [0 1 0; 1 0 1; 0 1 0];
sigy3 = [0 -1i 0; i 0 -1i; 0 i 0];
sigz3 = [1 0 0; 0 0 0; 0 0 -1];

Sx_1 = 1/sqrt(2)*sigx3;
Sx_2 = 0.5.*sigx;
Sx_22 = 0.5.*sigx;

Sy_1 = 1/sqrt(2)*sigy3;
Sy_2 = 0.5.*sigy;
Sy_22 = 0.5.*sigy;

Sz_1 = sigz3;
Sz_2 = 0.5.*sigz;
Sz_22 = 0.5.*sigz;

% limit = 100; %set iteration limit
D = 1; %setting dimensions
k = 0.0862875; %Boltzmann constant
J_T = 20.0;
J_J = 2.0;
minField = -0.29;
maxField = 0.29;
fieldStep = 5e-4;

% minTemp = 0.01;
% maxTemp = 1.51;
% tempStep = 0.5;
% temperature = int8((maxTemp-minTemp)/tempStep+1); % discretization of the temperature
% tempSeries = [minTemp:tempStep:maxTemp];

tempSeries = [0.01, 0.02, 0.05, 0.2, 1.0];
temperature = length(tempSeries); % Alternative temperature steps
J_1 = J_T*(k*tempSeries(1)); %coupling strength
J_T = [J_1./(k.*tempSeries)];
J_2 = J_J * J_1;

field = int8((maxField-minField)/fieldStep); % discretization of the field
%vector = zeros(int8(field*temperature),3);
% map = zeros(field, temperature);

delta = 0.0001;
iterator = 1;
newAvSx_1 = 0;
newAvSx_2 = 0;
newAvSx_22 = 0;
avmx_up = zeros(temperature, field);
avmx_down = zeros(temperature, field);

newAvSy_1 = 0;
newAvSy_2 = 0;
newAvSy_22 = 0;
avmy_up = zeros(temperature, field);
avmy_down = zeros(temperature, field);

newAvSz_1 = 0;
newAvSz_2 = 0;
newAvSz_22 = 0;
avmz_up = zeros(temperature, field);
avmz_down = zeros(temperature, field);

ti = 1;
for temp = tempSeries
    
    beta = 1/(k*temp);
    
    avSx_1 = 0; %initial guess of average Sx(1)
    avSx_2 = 0; %initial guess of average Sx(1/2)
    avSx_22 = 0; %initial guess of average Sx(1/2)_2
    
    avSy_1 = 0; %initial guess of average Sy(1)
    avSy_2 = 0.5; %initial guess of average Sy(1/2)
    avSy_22 = -0.5; %initial guess of average Sy(1/2)_2
    
    avSz_1 = 1; %initial guess of average Sz(1)
    avSz_2 = 0.5; %initial guess of average Sz(1/2)
    avSz_22 = 0.5; %initial guess of average Sz(1/2)_2
    
%     avSx = avSx_1 + avSx_2 + avSx_22;
%     avSy = avSy_1 + avSy_2 + avSy_22;
%     avSz = avSz_1 + avSz_2 + avSz_22;
    
    hi = 1;
    for hL = minField:fieldStep:maxField
        while 1
            %             Hx = J_1*avSx_2*kron(kron(Sx_1,eye(2)),eye(2)) + J_2*avSx_2*kron(kron(Sx_2,eye(2)),eye(3));
            %             Hy = J_1*avSy_2*kron(kron(Sy_1,eye(2)),eye(2)) + J_2*avSy_2*kron(kron(Sy_2,eye(2)),eye(3));
            %             Hz = J_1*avSz_2*kron(kron(Sz_1,eye(2)),eye(2)) + J_2*avSz_2*kron(kron(Sz_2,eye(2)),eye(3));
            
            Hx_12 = J_1*avSx_2*kron(kron(Sx_1,eye(2)),eye(2));
            Hx_12_2 = J_1*avSx_22*kron(kron(Sx_1,eye(2)),eye(2));
            Hx_22 = J_2*avSx_22*kron(kron(Sx_2,eye(2)),eye(3));
            Hx_22_2 = J_2*avSx_2*kron(kron(Sx_22,eye(2)),eye(3));
            
            Hy_12 = J_1*avSy_2*kron(kron(Sy_1,eye(2)),eye(2));
            Hy_12_2 = J_1*avSy_22*kron(kron(Sy_1,eye(2)),eye(2));
            Hy_22 = J_2*avSy_22*kron(kron(Sy_2,eye(2)),eye(3));
            Hy_22_2 = J_2*avSy_2*kron(kron(Sy_22,eye(2)),eye(3));
            
            Hz_12 = J_1*avSz_2*kron(kron(Sz_1,eye(2)),eye(2));
            Hz_12_2 = J_1*avSz_22*kron(kron(Sz_1,eye(2)),eye(2));
            Hz_22 = J_2*avSz_22*kron(kron(Sz_2,eye(2)),eye(3));
            Hz_22_2 = J_2*avSz_2*kron(kron(Sz_22,eye(2)),eye(3));
            
            % Further assemble the hamiltonian and apply 0.5 factor to
            % avoid double counting
            Hx = Hx_12 + Hx_12_2 + 0.5*Hx_22 + 0.5* Hx_22_2;
            Hy = Hy_12 + Hy_12_2 + 0.5*Hy_22 + 0.5* Hy_22_2;
            Hz = Hz_12 + Hz_12_2 + 0.5*Hz_22 + 0.5* Hz_22_2;
            
            Hamlt = Hx + Hy + Hz - hL*kron(eye(4),Sz_1) - hL*kron(eye(6),Sz_2) - hL*kron(eye(6),Sz_22); %with external field
            %              Hamlt = Hx + Hy + Hz; %without external field
            [v, E] = eig(Hamlt);
            Z = trace(exp(-beta.*E));
            newAvSx_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sx_1)*v)/Z;
            newAvSx_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sx_2)*v)/Z;
            newAvSx_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sx_22)*v)/Z;
            diffx_1 = newAvSx_1 - avSx_1;
            diffx_2 = newAvSx_2 - avSx_2;
            diffx_22 = newAvSx_22 - avSx_22;
            
            newAvSy_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sy_1)*v)/Z;
            newAvSy_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sy_2)*v)/Z;
            newAvSy_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sy_22)*v)/Z;
            diffy_1 = newAvSy_1 - avSy_1;
            diffy_2 = newAvSy_2 - avSy_2;
            diffy_22 = newAvSy_22 - avSy_22;
            
            newAvSz_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sz_1)*v)/Z;
            newAvSz_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sz_2)*v)/Z;
            newAvSz_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sz_22)*v)/Z;
            diffz_1 = newAvSz_1 - avSz_1;
            diffz_2 = newAvSz_2 - avSz_2;
            diffz_22 = newAvSz_22 - avSz_22;
            
            if abs(diffx_1) < delta && abs(diffx_2) < delta && abs(diffy_1) < delta && abs(diffy_2) < delta && abs(diffz_1) < delta && abs(diffz_2) < delta && abs(diffz_22) < delta
                break
            else
                avSx_1 = newAvSx_1;
                avSx_2 = newAvSx_2;
                avSx_22 = newAvSx_22;
                
                avSy_1 = newAvSy_1;
                avSy_2 = newAvSy_2;
                avSy_22 = newAvSy_22;
                
                avSz_1 = newAvSz_1;
                avSz_2 = newAvSz_2;
                avSz_22 = newAvSz_22;
            end
        end
        avmx_up(ti,hi) = avSx_1 + avSx_2 + avSx_22;
        avmy_up(ti,hi) = avSy_1 + avSy_2 + avSy_22;
        avmz_up(ti,hi) = avSz_1 + avSz_2 + avSz_22;
        
        hi = hi + 1;
    end
    
    avSx_1 = 0; %initial guess of average Sx(1)
    avSx_2 = 0; %initial guess of average Sx(1/2)
    
    avSy_1 = 0; %initial guess of average Sy(1)
    avSy_2 = 0; %initial guess of average Sy(1/2)
    
    avSz_1 = 0; %initial guess of average Sz(1)
    avSz_2 = 0; %initial guess of average Sz(1/2)
%     hi = 1;
%     for hL = maxField:-fieldStep:minField
%         while 1
%             %             Hx = J_1*avSx_2*kron(kron(Sx_1,eye(2)),eye(2)) + J_2*avSx_2*kron(kron(Sx_2,eye(2)),eye(3));
%             %             Hy = J_1*avSy_2*kron(kron(Sy_1,eye(2)),eye(2)) + J_2*avSy_2*kron(kron(Sy_2,eye(2)),eye(3));
%             %             Hz = J_1*avSz_2*kron(kron(Sz_1,eye(2)),eye(2)) + J_2*avSz_2*kron(kron(Sz_2,eye(2)),eye(3));
%             
%             Hx_12 = J_1*avSx_2*kron(kron(Sx_1,eye(2)),eye(2));
%             Hx_12_2 = J_1*avSx_22*kron(kron(Sx_1,eye(2)),eye(2));
%             Hx_22 = J_2*avSx_22*kron(kron(Sx_2,eye(2)),eye(3));
%             Hx_22_2 = J_2*avSx_2*kron(kron(Sx_22,eye(2)),eye(3));
%             
%             Hy_12 = J_1*avSy_2*kron(kron(Sy_1,eye(2)),eye(2));
%             Hy_12_2 = J_1*avSy_22*kron(kron(Sy_1,eye(2)),eye(2));
%             Hy_22 = J_2*avSy_22*kron(kron(Sy_2,eye(2)),eye(3));
%             Hy_22_2 = J_2*avSy_2*kron(kron(Sy_22,eye(2)),eye(3));
%             
%             Hz_12 = J_1*avSz_2*kron(kron(Sz_1,eye(2)),eye(2));
%             Hz_12_2 = J_1*avSz_22*kron(kron(Sz_1,eye(2)),eye(2));
%             Hz_22 = J_2*avSz_22*kron(kron(Sz_2,eye(2)),eye(3));
%             Hz_22_2 = J_2*avSz_2*kron(kron(Sz_22,eye(2)),eye(3));
%             
%             % Further assemble the hamiltonian and apply 0.5 factor to
%             % avoid double counting
%             Hx = Hx_12 + Hx_12_2 + 0.5*Hx_22 + 0.5* Hx_22_2;
%             Hy = Hy_12 + Hy_12_2 + 0.5*Hy_22 + 0.5* Hy_22_2;
%             Hz = Hz_12 + Hz_12_2 + 0.5*Hz_22 + 0.5* Hz_22_2;
%             
%             Hamlt = Hx + Hy + Hz - hL*kron(eye(4),Sz_1) - hL*kron(eye(6),Sz_2) - hL*kron(eye(6),Sz_22); %with external field
%             %              Hamlt = Hx + Hy + Hz; %without external field
%             [v, E] = eig(Hamlt,'nobalance');
%             Z = trace(exp(-beta.*E));
%             newAvSx_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sx_1)*v)/Z;
%             newAvSx_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sx_2)*v)/Z;
%             newAvSx_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sx_22)*v)/Z;
%             diffx_1 = newAvSx_1 - avSx_1;
%             diffx_2 = newAvSx_2 - avSx_2;
%             diffx_22 = newAvSx_22 - avSx_22;
%             
%             newAvSy_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sy_1)*v)/Z;
%             newAvSy_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sy_2)*v)/Z;
%             newAvSy_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sy_22)*v)/Z;
%             diffy_1 = newAvSy_1 - avSy_1;
%             diffy_2 = newAvSy_2 - avSy_2;
%             diffy_22 = newAvSy_22 - avSy_22;
%             
%             newAvSz_1 = diag(exp(-beta*E))'*diag(v'*kron(eye(4),Sz_1)*v)/Z;
%             newAvSz_2 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sz_2)*v)/Z;
%             newAvSz_22 = diag(exp(-beta*E))'*diag(v'*kron(eye(6),Sz_22)*v)/Z;
%             diffz_1 = newAvSz_1 - avSz_1;
%             diffz_2 = newAvSz_2 - avSz_2;
%             diffz_22 = newAvSz_22 - avSz_22;
%             
%             if abs(diffx_1) < delta && abs(diffx_2) < delta && abs(diffy_1) < delta && abs(diffy_2) < delta && abs(diffz_1) < delta && abs(diffz_2) < delta && abs(diffz_22) < delta
%                 break
%             else
%                 avSx_1 = newAvSx_1;
%                 avSx_2 = newAvSx_2;
%                 avSx_22 = newAvSx_22;
%                 
%                 avSy_1 = newAvSy_1;
%                 avSy_2 = newAvSy_2;
%                 avSy_22 = newAvSy_22;
%                 
%                 avSz_1 = newAvSz_1;
%                 avSz_2 = newAvSz_2;
%                 avSz_22 = newAvSz_22;
%             end
%         end
%         avmx_down(ti,hi) = avSx_1 + avSx_2 + avSx_22;
%         avmy_down(ti,hi) = avSy_1 + avSy_2 + avSy_22;
%         avmz_down(ti,hi) = avSz_1 + avSz_2 + avSz_22;
%         
%         %         vector(iterator,:) = [temp, hL, newAvSz_1 + 2 * newAvSz_2];
%         %         map(floor(iterator/(temperature+1))+1,rem(iterator,temperature+1)+1) = avm;
%         iterator = iterator + 1;
%         %         if mod(iterator, 10) == 0
%         %            sprintf('Current iteration: %d, and Current temperature: %d', iterator, temp);
%         %         end
%         hi = hi + 1;
%     end
    ti = ti + 1;
end

figure
hold on
grid on
axis([1.2*minField 1.2*maxField -2 2]);
xlabel('External magnetic filed (H)');
ylabel('Magnetization per unit cell (<m>)');
for i = 1:temperature
    plot(minField:fieldStep:maxField, avmx_up(i,:),'o-','MarkerSize',5);
%     plot(maxField:-fieldStep:minField, avmx_down(i,:),'s-');
end
hold off

figure
hold on
grid on
axis([1.2*minField 1.2*maxField -2 2]);
xlabel('External magnetic filed (H)');
ylabel('Magnetization per unit cell (<m>)');
for i = 1:temperature
    plot(minField:fieldStep:maxField, avmy_up(i,:),'o-','MarkerSize',5);
%     plot(maxField:-fieldStep:minField, avmy_down(i,:),'s-','MarkerSize',5);
end
hold off

figure
hold on
grid on
axis([1.2*minField 1.2*maxField -2 2]);
xlabel('External magnetic filed (H)');
ylabel('Magnetization per unit cell (<m>)');
for i = 1:temperature
    plot(minField:fieldStep:maxField, avmz_up(i,:),'o-','MarkerSize',5);
    %plot(maxField:-fieldStep:minField, avmz_down(i,:),'s-','MarkerSize',5);
end

lgd = strcat('T = ', string(num2cell((tempSeries))), ' K, J_1/k_BT = ', string(num2cell(J_T)), ', J_2/J_1 = ', num2str(J_J));
legend(lgd);
hold off
% plot(vector(:,2),vector(:,3));
% scatter3(vector(:,1),vector(:,2),vector(:,3),15,vector(:,3),'filled');
% view([0 90]);