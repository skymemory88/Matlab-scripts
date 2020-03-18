%% Caluculate the timedependence of the magnetic order in CoCl2 2D2O
% Starting from the ferrimagnetic phase long range antiferromagnetic order will 
% By Ursula Bengaard Hansen, 23/07/2019



%%
clear
close all
home;
tic
rng shuffle
%%
grey    = [0.6 0.6 0.6];
colAF   = [0.3922 0.5843 0.9294];
colFI   = [0.5412 0.1686 0.8863];
colF    = [0.7804    0.0824    0.5216];

colmap = colormap('jet');
%%
factor = 1;
plot_key = [1 1];
if plot_key(1) 
    %Np = [1e3 5e3 1e4 1.5e4 5e5 1e6 1.5e6 5e6].*factor.^2;
    Np = [1e3 5e3 1.5e4 5e5 1e6 1.5e6 5e6].*factor.^2;
else
    Np  = [1 500 1000 2500:2500:2e4 3e4:1e4:1e5 2e5:1e5:5e6].*factor.^2;
end

L   = 60.*factor;        % Laticce dimentions has to be even

fprintf('Lattice dimensions : %2.2f x %2.2f \n\n',L,L);

h   = 0;         % gmuBH/J1z
T   = 1e-4;      % KbT/J1z

fprintf('Field : %2.2e \n\n',h);
fprintf('Temperature : %2.2e \n\n',T);

J1z = 0.4029;    % meV
J1  = 1;         % J1z/J1z
J2  = 0.2154;    % J2z/J1z

J = [J1 J2 0];
spin = 0.5;

N   = Np(end);
dim = [L L];
Fdim = [2*L 2*L];
diaFSNp = zeros(length(Np),Fdim(2));

% apodization window 
w=window2(L,L,@hann);
%%
% initial Ferrimagnetic spin configuration
Si       = repmat([-1 1 1; 1 1 -1;1 -1 1 ],floor(dim(1)/3),floor(dim(2)/3));
[Ei,Si] = calculate_energy(Si,h,J,spin,dim);
E1 = Ei; S1 = Si;
FS1 = fft2(S1(1:end-1,1:end-1).*w,Fdim(1),Fdim(2));
FS1 = abs(FS1);
diaFS1 = diag(FS1)./(Fdim(1)*Fdim(2));
FS12 = log(FS1);

% final AF configuration
SAF = repmat([-1 1; 1 -1],floor(dim(1)/2),floor(dim(2)/2));
[EAF,SAF] = calculate_energy(SAF,h,J,spin,dim);
FSAF = fft2(SAF(1:end-1,1:end-1).*w,Fdim(1),Fdim(2));
FSAF = abs(FSAF);
diaFSAF = diag(FSAF)./(Fdim(1)*Fdim(2));
FSAF2 = log(FSAF);

%% Calculate time dependence using Monte Carlo rutine

m = 1;
prog = 10;
for n1=1:N
    % flip one spin
    i = randi(dim(1)); j = randi(dim(2));
    Sn      = Si(1:end-1,1:end-1);
    Sn(i,j) = -Sn(i,j);
    
    % Calculate the energy of the new configuration
    [En,Sn] = calculate_energy(Sn,h,J,spin,dim);
    
    dE = En - Ei;
    
    % Monte Carlo code
    if rand <= min(exp(-dE/T),1) % Accept the sample with prob = min(exp(-dE/T),1)
        Si = Sn;
        Ei = En;
    end
    
    % Plot spin state at time Np
    if ismember(n1,Np)
        
        FSn = fft2(Si(1:end-1,1:end-1).*w,Fdim(1),Fdim(2));
        FSn = abs(FSn);
        FSn2 = log(FSn);
        
        if plot_key(1)
            figure
            imagesc(FSn2);
            caxis([0 8])
            colormap('parula')
%             export_fig(['FS_L' num2str(L) '_N' num2str(n1) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf','-OpenGl')
            
            figure
            imagesc(Si(1:end-1,1:end-1));
            colormap('gray')
%             export_fig(['S_L' num2str(L) '_N' num2str(n1) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf','-OpenGl')
        end
        
        diaFSNp(m,:) = diag(FSn)./(Fdim(1)*Fdim(2));
        m=m+1;
    end
    
    % write out progress to terminal
    if mod(n1,N/10) == 0
        fprintf('%1.0f %% \t',prog);
        prog = prog +10;
    end
end

 fprintf('\n\n');
 toc
Sf = Si(1:end-1,1:end-1);
Ef = Ei;
dE = Ef-E1;

%% Final configuration
FSf = fft2(Sf.*w,Fdim(1),Fdim(2));
FSf = abs(FSf);
FSf2 = log(FSf);

if plot_key(2)
figure
imagesc(FSf2);
caxis([0 8])
% export_fig(['FS_L' num2str(L) '_N' num2str(N) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf','-OpenGl')

figure
imagesc(Sf) 
colormap('gray')
% export_fig(['S_L' num2str(L) '_N' num2str(N) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf','-OpenGl')
end

%% Q integrated
col = makeColorMap([1 0 0],[0 0 1],length(Np));
figure
hold on
diaX = 2.*(0:1/Fdim(2):(Fdim(2)-1)/Fdim(2));

plot(diaX,diaFS1,'-k')
plot(diaX,diaFSAF,'-k')
for i = 1:length(Np)
    plot(diaX,diaFSNp(i,:)','-','Color',col(i,:))
end
xlabel('[H 0 0]')
hold off
% export_fig(['I_H00_L' num2str(L) '_N' num2str(N) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf')

%% Plot time dependence
iAF = floor(Fdim(1)/2+1);
iFI = floor(Fdim(1)/3+1);

n1 = (log(diaFSNp(end,iAF))-log(diaFSNp(1,iAF)))/(log(Np(end))-log(Np(1)));
A1 = exp(log(diaFSNp(1,iAF))-n1*log(Np(1)));

n2 = 0.5;
A2 = exp(log(diaFSNp(1,iAF))-n2*log(Np(1)));

figure
hold on
plot([0 Np(end)],[diaFS1(iFI) diaFS1(iFI)],'-.','Color',colFI)
plot([0 Np(end)],[diaFSAF(iAF) diaFSAF(iAF)],'-.','Color',colAF)

plot(Np(:),A1*Np(:).^(n1),':k')
plot(Np(:),A2*Np(:).^n2,'--k')

hFI=plot(Np(1:end),diaFSNp(1:end,iFI),'o-','Color',colFI);
hAF=plot(Np(1:end),diaFSNp(1:end,iAF),'o-','Color',colAF);

ylim([0 0.1])
xlabel('Number of steps')
legend([hAF hFI],'AF','FI','location','northwest')
% export_fig(['I_N_L' num2str(L) '_N' num2str(N) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf')

%% logscale
figure
ax_log = gca;
set(ax_log,'YScale','log');
set(ax_log,'XScale','log');
set(ax_log,'YMinorTick','on','XMinorTick','on');

hold on
plot([1 Np(end)],[diaFS1(iFI) diaFS1(iFI)],'-.','Color',colFI)
plot([1 Np(end)],[diaFSAF(iAF) diaFSAF(iAF)],'-.','Color',colAF)

plot([1 Np(end)],[A1*1^n1 A1*Np(end)^n1],':k')
plot([1 Np(end)],[A2*1^n2 A2*Np(end)^n2],'--k')

hFI=plot(Np(1:end),diaFSNp(1:end,iFI),'o-','Color',colFI);
hAF=plot(Np(1:end),diaFSNp(1:end,iAF),'o-','Color',colAF);

xlabel('Number of steps')
xlim([1 Np(end)])
% ylim([1e-4 1])
legend([hAF hFI],'AF','FI','location','southeast')

%%
figure
hold on
plot([0 Np(end)],[diaFS1(iFI) diaFS1(iFI)],'-.','Color',colFI)
plot([0 Np(end)],[diaFSAF(iAF) diaFSAF(iAF)],'-.','Color',colAF)
hFI=plot(Np(1:end),diaFSNp(1:end,iFI),'o-','Color',colFI);
hAF=plot(Np(1:end),diaFSNp(1:end,iAF),'o-','Color',colAF);
ylim([0 0.1])
xlabel('Number of steps')
xlim([0 2.5e6].*factor.*factor )
legend([hAF hFI],'AF','FI','location','northwest')
% export_fig(['I_Nzoom_L' num2str(L) '_N' num2str(N) '_H' num2str(h*10) '_T' num2str(T*1e5) 'E-5_'],'-pdf')

%%
diaTotal = [diaFS1'; diaFSAF'; diaFSNp];
% save 'output.mat' diaTotal Np L h T;

toc
