% B = [-57.9   0.309   3.51   0.00   0.000540   0.0631   0.0171]; % Ho -- Phys. Rev. B 92, 144422 (2015)
% B = [-60.0   0.350   3.60   0.00   0.000400   0.0700   0.0098]; % Ho -- Phys. Rev. B 75, 054426 (2007)
% B = [-60.0   0.350   3.60   0.00   0.000400   0.0655   0.0098]; % CEF_2007_mod
B = [-57.9   0.309   3.60   0.00   0.000540   0.0570   0.0130]; % CEF_2015_mod
% B = [-54.7   0.262   2.91   0.00   0.000400   0.0511   0.0138]; % Ho -- PRB 2015 low limit 
% B = [-61.1   0.356   4.11   0.00   0.000680   0.0751   0.0204]; % Ho -- PRB 2015 high limit 
% B = [-57.9   0.309   3.60   0.00   0.000540   0.0565   0.0165]; % Ho -- SC239 (Jz/1) {phi=0, nZ=0, demag=1}
% B = [-57.9   0.309   3.51   0.00   0.000540   0.0560   0.0148]; % Ho -- SC239 (Jz/1) {phi=0.5, nZ=0, demag=0}
temp = 4; % tempertaure
rot = 11; % crystal field rotation relative to the crystallographic axes

% Deduce the crystal-field Hamiltonian
J=8;
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-(J-1:-1:-J)).*(J+1+(J-1:-1:-J))),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

% Diagonalise the crystal-field Hamiltonian
Hcf = cf(J, B/1e3, rot); % crystal field energy only
[wav, ee] = eig(Hcf);
ee = real(diag(ee));            
En = ee - min(ee); % normalize to the ground state

% opt.ion = 'Ho';
% opt.eUnit = 'meV';
% opt.plot = false;
% opt.ScanMode = 'temp';
% opt.hyp = 0.0;
% [res, ~] = EngyLevels(0,0,0,17,opt);
% En = res.En{:};
% En = squeeze(En(1,:));
% wav = res.wav(:,:);

% Calculate the transition matrix elements
Z = sum(exp(-11.6*En/temp)); % Define the partition function
ll = 1;
for ni = 1:8 % initial state
    for nf = 1:8 % final state
        if ni ~= nf && round(1e6*En(nf)) ~= round(1e6*En(ni)) && En(nf) > En(ni)
            aJx = wav(:,nf)'*Jx*wav(:,ni);
            aJy = wav(:,nf)'*Jy*wav(:,ni);
            aJz = wav(:,nf)'*Jz*wav(:,ni);

            nm(ll,:) = [ni nf];                                     % Matrix element considered
            dE(ll) = En(nf)-En(ni);                                    % Energy difference
            S(ll) = ( dot(aJx,aJx) + dot(aJy,aJy) + dot(aJz,aJz) )*...
                exp(-11.6*En(ni)/temp)/Z; % S(q,w), intensity
            ll = ll+1;
        end
    end
end
    

[dE, ia, ic] = unique(round(1e6*dE)); % round to keV
dE = dE/1e6; % [meV]
for n = 1:max(ic)        
    Sqw(n) = sum(S(ic == n));           % Sum intensities for a given transition
    multi(n) = sum(ic == n);
end  

% Number of crystal field levels   
ngau = length(dE);

amp(1:3:3*ngau) = Sqw;        % define amplitudes of peaks
amp(2:3:3*ngau) = dE;         % define crystal field energies
amp(3:3:3*ngau) = 0.2;        % define widths of peaks
amp(3*ngau+1) = 0;

vp(1:3:3*ngau) = 0;         % define amplitudes of peaks
vp(2:3:3*ngau) = 0;         % define crystal field energies
vp(3:3:3*ngau) = 0;         % define widths of peaks
vp(3*ngau+1) = 0;

Es = 0:0.01:8.7; %[meV]

ngss = (length(amp)-1)/3;
sc = zeros(size(Es));
for ii = 1:ngss
    off = 3*(ii-1);
    sc = sc + amp(1+off)*exp(-0.5*(((Es-amp(2+off))/amp(3+off)).^2)/sqrt(2*pi)/amp(3+off));
end
sc = sc + amp(4+off);
sc = sc./10; % arbituray unit change

%% Plot results ===========================================================
hfig = figure(3234);
% clf
plot(Es,sc,'b')
hold on
title([num2str(temp) ' K'])
xlabel('E (meV)')
ylabel('S(Q,\omega) (arb.)')

% % hl(1) = line([0 0],[-0.1 0.7],'color','k');
% for n = 1:ngss
%     hl(n) = line([dE(n) dE(n)],[0 12],'color','k');
% end
% xlim([-0.7 8.7])
