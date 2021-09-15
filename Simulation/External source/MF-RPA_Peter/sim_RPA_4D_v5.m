function varargout = sim_RPA_4D_v5(sample_in)

% If there is an input, use that to overwrite the default parameters
% This code is the first to use J.Jensen's RPA for binary alloy systems for
% calculating S(q,w)

%% Calling following functions (not exhaustive list):
%     'LiErxHoyY1_x_yF4_mod'
%     'chi_qw_VC'
%     'remf_mod'
%     'dipole_direct'
%     'cf_Ho'
%     'cf_Er'
%     'exchange'
%     'ellipsoid_demagn'

%% Version build: 14.03.2018 by PB

% Since 08.06.2017:
% - Changed the way chi_qw_VC calculates chi00 ... chi11 (eliminated b*inv(M)
% operation)
% - Added a check if one of the domain populations is 0
% - Corrected the (delta - q.q) dipolar factor to be in A-1 not in rlu
% 29.10.2018 / v5: 
% - Save chi1,chi2,chitot tensor as a (T, H, 3, 3, Q, w)
% - included h4 anisotropy as h4(Jx^4 + Jy^4)

%% Setup parameters for MF calculation -------------------------------------
    
% Remove all global variables
clearvars -global

global strategies; % global convergence strategies switches
strategies.powerlaw = false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator = 0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping = 0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit = true; % turn on fitting to an exponential.
strategies.expfit_period = 30; % period between exponential fit.
strategies.expfit_deltaN = 5; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.
% -------------------------------------------------------------------------

%% Definitions for the simulations
% Field or range of fields
sample.hvec = linspace(0,0.7,201);
sample.fields=[zeros(size(sample.hvec))
               zeros(size(sample.hvec))
               sample.hvec];

% Temperature (or range of temperatures)
sample.temp = 0.1;

sample.Erbium = 1.0; % proportion erbium
sample.Holmium = 1-sample.Erbium; % proportion holmium

% Super-exchange interactions
sample.exEr = 0; % magnetic exchange in J*S_iS_j
sample.exHo = -0.0001; % 0.1 microeV from Henrik and Jens PRB

% Hyperfine coupling parameters
sample.A_Er = 0.00043412;   % hyperfine coupling energy in meV for Erbium (Phys. Rev. B 2, 2298 - 2301 (1970) )
sample.A_Ho = 0.003361;     % hyperfine coupling energy in meV for Holmium (from Phys. Rev. B 75, 054426 (2007) ).

% Renormalization factors
sample.renorm_Ho = [1,1,0.785];
sample.renorm_Er = [1,1,1];

% Proportion of RE ions carrying a nuclear spin
sample.ErHyp = 0.23; % proportion, whithin the erbium percentage, of isotops carrying a nuclear spin (23% is default).
sample.HoHyp = 1;

% Include or not the effect of the demagnetisation factor
sample.demagn = true;

% Sample shape
sample.alpha = 0; % shape is a sphere (1), needle (0), or disk (inf)

%J and gLande
sample.J_Er = 15/2;
sample.gLande_Er = 1.2;
sample.I_Er = 7/2;
sample.J_Ho = 8;
sample.gLande_Ho = 1.25;
sample.I_Ho = 7/2;

% S(q,omega) range:
% sample.omega = linspace(-0.01,0.25,201);     % meV E-cut
sample.epsilon = 0.002;                       % meV damping term
% sample.epsilon = 1.5e-4;
sample.pop = [0.5 0.5];                       % domain population

% Define where in reciprocal space to simulate
sample.qvec = 0; %linspace(-1,0,7);
sample.h = 0;
sample.k = 0;
sample.l = 0;

% Define lattice parameters:
% sample.abc = [5.175 0 0
%               0 5.175 0
%               0 0 10.75]; % Ho
sample.abc = [5.162 0 0
              0 5.162 0
              0 0 10.70]; % Er
        
% Unit cell volume:
sample.vol = sum(sample.abc(1,:).*cross(sample.abc(2,:),sample.abc(3,:)));
% Reciproca l lattice unit vectors
sample.ABC = [2*pi*cross(sample.abc(2,:),sample.abc(3,:))/sample.vol
              2*pi*cross(sample.abc(3,:),sample.abc(1,:))/sample.vol
              2*pi*cross(sample.abc(1,:),sample.abc(2,:))/sample.vol];                   

%% For testing and bug-proofing purposes
sample.mixhyp = true;                       % calculate as a binary compound with both HF and non-HF ions
sample.hyp_nohyp = false;                   % replace chi0 in VC calculation to be the same
sample.nohyp_hyp = false;                   % Replace chi(no HF) with chi(HF)
sample.verdate = '2018-10-29';              % date of the version of the RPA program used
sample.ver = 'v5';                          % version number
sample.rmHam = true;                        % the final "sample" structure does not contain the calculated Ham
sample.rmint = false;                       % calculate only 1/(1 - Jx) in the dynamic susceptibility (remove intensities)
% sample.h4 = 0;                        % introduce h4 anisotropy (meV) due to order by disorder, h4(Jx^4 + Jy^4)
%% Set input parameters defined in sample_in structure
% (sim_RPA_4D_v5 can be called by other functions)

if nargin > 0
    fnames = fieldnames(sample_in);
    for n = 1:length(fnames)
        sample.(fnames{n}) = sample_in.(fnames{n});
    end
end

%% Tidy up input
q = makeQ(sample.h,sample.k,sample.l);

% Define sizes of arrays
Nq = size(q,1);
Ne = length(sample.omega);
Nh = length(sample.hvec);
Nt = length(sample.temp);

sample.S1 = zeros(Nt,Nh,Nq,Ne);            % Domain 1
sample.S2 = zeros(Nt,Nh,Nq,Ne);            % Domain 2
sample.Stot = zeros(Nt,Nh,Nq,Ne);          % Domain averaged
sample.moment0 = zeros(Nt,Nh,3);       % keep the moment of the 1st ion
sample.chi1 = zeros(Nt,Nh,3,3,Nq,Ne);            % Domain 1
sample.chi2 = zeros(Nt,Nh,3,3,Nq,Ne);            % Domain 2
sample.chitot = zeros(Nt,Nh,3,3,Nq,Ne);            % Domain averaged

%% Display input parameters:

warning off

Hx = sample.fields(1,:);
Hy = sample.fields(2,:);
Hz = sample.fields(3,:);

if sample.Erbium == 1
    sample.ionname = 'Er';
end
if sample.Holmium == 1
    sample.ionname = 'Ho';
end

disp('===========================================================================================')
disp(['Simulation parameters for LiEr(' num2str(sample.Erbium) ')Ho(' num2str(sample.Holmium) ')F4'])
if length(sample.temp) > 1
    disp([' T  (K) = ' num2str([min(sample.temp) max(sample.temp) sample.temp(2)-sample.temp(1) ], '%3.2f-%3.2f, dT = %3.2f')])
else
    disp([' T  (K) = ' num2str(sample.temp, '%3.2f')])
end
if length(Hx) > 1
    disp([' Hx (T) = ' num2str([min(Hx) max(Hx) Hx(2)-Hx(1) ], '%3.2f-%3.2f, dHx = %3.2f')])
else
    disp([' Hx (T) = ' num2str(Hx, '%3.2f')])
end
if length(Hy) > 1
    disp([' Hy (T) = ' num2str([min(Hy) max(Hy) Hy(2)-Hy(1) ], '%3.2f-%3.2f, dHx = %3.2f')])
else
    disp([' Hy (T) = ' num2str(Hy, '%3.2f')])
end
if length(Hz) > 1
    disp([' Hz (T) = ' num2str([min(Hz) max(Hz) Hz(2)-Hz(1) ], '%3.2f-%3.2f, dHx = %3.2f')])
else
    disp([' Hz (T) = ' num2str(Hz, '%3.2f')])
end
disp(num2str(sample.(['J_' sample.ionname]), ' J = %3.2f'))
% disp(num2str(ion.B(nion,:)*1e3,' B_CEF (meV) = [%3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f]'))
disp(num2str([sample.([sample.ionname 'Hyp']) sample.(['I_' sample.ionname]) sample.(['A_' sample.ionname])*1e3 sample.(['ex' sample.ionname])*1e6],'f = %3.2f, I = %3.3f, A (ueV) = %3.3f, J_ex (ueV) = %3.3f'))
% disp('Initial moment configuration: ')
% disp(num2str(ion.mom(:,:,nion),'%3.2f, %3.2f, %3.2f\n %3.2f, %3.2f, %3.2f\n %3.2f, %3.2f, %3.2f'))
%disp(num2str(sample.mixhyp,'If f is non-zero, Sqw = fSqw(hyp) + (1-f)Sqw(nohyp)? : %d'))
disp(num2str([sample.demagn sample.alpha],'Demagnetisation included: %d, Sample shape: %d'))
disp(num2str(sample.(['renorm_' sample.ionname]),'Renormalisation factor: %3.2f, %3.2f, %3.2f'))
disp(num2str([sample.rmHam sample.rmint sample.h4*1000],'rm Ham: %d, rm Intensity: %d, h4 (ueV) = %3.3f'))
disp('--------------------------------------------------------------')
if length(sample.h) > 1
    disp([' h = ' num2str([min(sample.h) max(sample.h) sample.h(2)-sample.h(1) ], '%3.2f->%3.2f, dh = %3.2f')])
else
    disp([' h = ' num2str(sample.h, '%3.2f')])
end
if length(sample.k) > 1
    disp([' k = ' num2str([min(sample.k) max(sample.k) sample.k(2)-sample.k(1) ], '%3.2f->%3.2f, dk = %3.2f')])
else
    disp([' k = ' num2str(sample.k, '%3.2f')])
end
if length(sample.l) > 1
    disp([' l = ' num2str([min(sample.l) max(sample.l) sample.l(2)-sample.l(1) ], '%3.2f->%3.2f, dl = %3.2f')])
else
    disp([' l = ' num2str(sample.l, '%3.2f')])
end
if Ne > 1
    disp([' w = ' num2str([min(sample.omega) max(sample.omega) sample.omega(2)-sample.omega(1) ], '%3.2f->%3.2f, dw = %3.6f meV')])
else
    disp([' w = ' num2str(sample.omega, '%3.2f meV')])
end
disp(num2str(sample.epsilon,'Damping = %3.7f meV'))
disp(num2str(sample.pop*100,'Domain population = %3.2f-%3.2f'))
disp(num2str([Nt Nh Nq Ne],'N(temp) = %d, N(field) = %d, N(Q) = %d, N(E) = %d'))
disp('=ver=29.10.2018============================================================================')

warning on


%% Calculate MF
sample = LiErxHoyY1_x_yF4_mod(sample);

%% Calculate MF-RPA
disp('... Running RPA')

if sample.pop(1) < 1e-4
    disp(['---> Skip calculation for domain 1, pop = ' num2str(sample.pop(1))])
end
if sample.pop(2) < 1e-4
    disp(['---> Skip calculation for domain 2, pop = ' num2str(sample.pop(2))])
end

tic
ll = 1;
for nh=1:Nh
    for nt=1:Nt
        
        chi0.hyp = zeros(3,3,4,Ne);
        chi0.nohyp = zeros(3,3,4,Ne);
        
        hh = sample.fields(:,nh)';
        t = sample.temp(nt);
        
        moment = sample.Jmom(:,:,nt,nh);
        if Nq == 1
            disp(sprintf(' %d of %d| %3.2f K [%3.2f %3.2f %3.2f] T | %3.2f %3.2f %3.2f, %3.2f %3.2f %3.2f, %3.2f %3.2f %3.2f, %3.2f %3.2f %3.2f',...
                ll,Nh*Nt,t, hh,reshape(moment',1,numel(moment))))
            ll = ll + 1;
        end
        sample.moment0(nt,nh,:) = moment(1,1:3);
        
        % For each ion in the unit cell, extract the mean-field Hamiltonian
        if sample.Erbium == 1
            for ionn = 1:4
                if sample.ErHyp ~= 0
                    Ham_hyp(ionn,:,:) = sample.Ham_hyp_Er_TH(ionn,nt,nh,:,:);
                end
                if sample.ErHyp == 0 || sample.mixhyp
                    Ham(ionn,:,:) = sample.Ham_Er_TH(ionn,nt,nh,:,:);
                end
            end
            
        end
        if sample.Holmium == 1
            for ionn = 1:4
                if sample.HoHyp ~= 0
                    Ham_hyp(ionn,:,:) = sample.Ham_hyp_Ho_TH(ionn,nt,nh,:,:);
                else
                    Ham(ionn,:,:) = sample.Ham_Ho_TH(ionn,nt,nh,:,:);
                end
            end
            
        end                
        
        %% Calculate single-ion chi0 from MF
        
        if sample.([sample.ionname 'Hyp']) ~= 0
            chi0.hyp = calc_chi0(sample,Ham_hyp,t,true);
            chi0.nohyp = zeros(size(chi0.hyp));
        else
            chi0.hyp = zeros(size(chi0.nohyp));
            chi0.nohyp = calc_chi0(sample,Ham,t,false);
        end

        
        %% Calculate S(q,w) for domain 1
        Skw1 = zeros(Nq,Ne);
        chi1 = zeros(3,3,Nq,Ne);
        
        if Nq > 1, ll = 1; end
        
        if sample.pop(1) > 1e-4           
            for nq = 1:Nq
                if Nq > 1
                    disp(sprintf(' %d of %d| %3.2f K [%3.2f %3.2f %3.2f] T | Q = (%3.2f %3.2f %3.2f)',...
                        ll,Nq*2,t, hh,q(nq,:)))
                    ll = ll + 1;
                end

                % calculate chi(q,w)
                chi1(:,:,nq,:) = chi_qw2(q(nq,:),chi0.hyp,sample); % added 2021.04.14 Yikai
%                 chi1(:,:,nq,:) = chi_qw(q(nq,:),chi0.hyp,sample); % Original code
%                 chi1(:,:,nq,:) = chi_qw_VC(q(nq,:),chi0,sample);

                Qvec = q(nq,:)*sample.ABC;     % convert into reciprocal lattice vectors
                
                for nw=1:length(sample.omega)
                    if ~sample.rmint
                        if t==0
                            Skw1(nq,nw)=1/pi* ...
                                sum(sum((eye(3)-(Qvec'*Qvec)/(Qvec*Qvec')).* ...
                                imag(chi1(:,:,nq,nw))));%-chim(:,:,nw,nq))));
                        else
                            Skw1(nq,nw)= 1/pi/(1-exp(-sample.omega(nw)*11.6/t))* ...
                                sum(sum((eye(3)-(Qvec'*Qvec)/(Qvec*Qvec')).* ...
                                imag(chi1(:,:,nq,nw))));%-chim(:,:,nw,nq))));
                        end
                    else
                        Skw1(nq,nw) = sum(sum(abs(imag(chi1(:,:,nq,nw)))));
                    end
                end

            end
        end
        
        sample.S1(nt,nh,:,:) = Skw1;
        sample.chi1(nt,nh,:,:,:,:) = chi1;
        
        %% Calculate S(q,w) for domain 2
        qq = [q(:,2) q(:,1) q(:,3)];
        
        % Remove q-points that are the same as calculated for domain #1
        indq = zeros(size(q,1),1);
        if sample.pop(1) > 1e-4 
            indq = sum((q - qq).^2,2)<1e-7;
        end
        
        % If input of q is a vector N x 3, calculate S(q,w) at each q point
        % without recalculating the entire single-ion Hamiltonian
        
        Skw2 = zeros(Nq,Ne);
        chi2 = zeros(3,3,Nq,Ne);
        
        if Nq > 1, ll = 1; end
        
        if sample.pop(2) > 1e-4            
            for nq = 1:Nq
                if Nq > 1
                    disp(sprintf(' %d of %d| %3.2f K [%3.2f %3.2f %3.2f] T | Q = (%3.2f %3.2f %3.2f)',...
                        Nq + ll,Nq*2,t, hh,q(nq,:)))
                    ll = ll + 1;
                end

                if indq(nq) == 0
                    % calculate chi(q,w)
                    chi2(:,:,nq,:) = chi_qw2(qq(nq,:),chi0.hyp,sample); % added 2021.04.14 Yikai
%                     chi2(:,:,nq,:) = chi_qw(qq(nq,:),chi0.hyp,sample); % Original code
%                     chi2(:,:,nq,:) = chi_qw_VC(qq(nq,:),chi0,sample);
                    
                    Qvec = qq(nq,:)*sample.ABC;     % convert into reciprocal lattice vectors
                    
                    for nw=1:length(sample.omega)
                        if ~sample.rmint
                            if sample.temp==0
                                Skw2(nq,nw)=1/pi* ...
                                    sum(sum((eye(3)-(Qvec'*Qvec)/(Qvec*Qvec')).* ...
                                    imag(chi2(:,:,nq,nw))));%-chim(:,:,nw,nq))));
                            else
                                Skw2(nq,nw)= 1/pi/(1-exp(-sample.omega(nw)*11.6/t))* ...
                                    sum(sum((eye(3)-(Qvec'*Qvec)/(Qvec*Qvec')).* ...
                                    imag(chi2(:,:,nq,nw))));%-chim(:,:,nw,nq))));
                            end
                        else
                            Skw2(nq,nw) = sum(sum(abs(imag(chi2(:,:,nq,nw)))));
                        end
                    end
                else
                    Skw2(nq,:) = Skw1(nq,:);
                end

            end        
        end
        
        sample.S2(nt,nh,:,:) = Skw2;
        sample.chi2(nt,nh,:,:,:,:) = chi2;
            
        % Total contribution to the scattering cross-section from both domains
        sample.Stot(nt,nh,:,:) = sample.pop(1)*Skw1 + sample.pop(2)*Skw2;
        sample.chitot(nt,nh,:,:,:,:) = sample.pop(1)*chi1 + sample.pop(2)*chi2;
        
    end
end
toc

% Convert meV->GHz (1 meV = 241.76811 GHz)
sample.meV2GHz = 241.76811;

if nargin == 0
    %% Make plots
    %
    % % 1) 2D plot: Field-Energy map of S(q,w)
    % figure(1003)
    % clf
    % hp = sanePColor(hh,omega,squeeze(sample.Stot)');
    % set(hp,'edgecolor','none')
    % xlabel('H_{\perp} (T)')
    % ylabel('omega (meV)')
    % colorbar
    % caxis([0 1e3])
    
    % 2) 2D plot:  Q-Energy map of S(q,w)
    %     figure(1005)
    %     clf
    %     hp = sanePColor(sample.qvec,sample.omega,squeeze(sample.Stot)');
    %     set(hp,'edgecolor','none')
    %     xlabel('(0,0,l) (rlu)')
    %     ylabel('omega (meV)')
    %     colorbar
    %     caxis([0 1e3])
    %% 3) 1D plot: omega-S(q,w)
    figure(1001)
    clf
    hp = plot(sample.omega,squeeze(sample.Stot)','-og');
    xlabel('omega (meV)')
    ylabel('S(q,w)')
end

%% EXIT -------------------------------------------------------------------
% I'm done, hurrah for that

if sample.rmHam
    sample = rmfield(sample,{'Ham_hyp_Er_TH','Ham_Er_TH','Ham_hyp_Ho_TH','Ham_Ho_TH',...
                                'Ham_hyp_Er','Ham_Er','Ham_hyp_Ho','Ham_Ho'});
end

if nargout > 0
    varargout{1} = sample;
end

end

function chi0 = calc_chi0(sample,Ham,t,hyper)

% Function to calculate chi0 for input into RPA, inc and exc HF

J = sample.(['J_' sample.ionname]);
I = sample.(['I_' sample.ionname]);

% Initiate J operators
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

if hyper == true
    Iid=eye(2*I+1);
%     Iz=diag(I:-1:-I);
%     hIz=kron(eye(2*J+1),Iz);
%     Ip=diag(sqrt((I-[(I-1):-1:-I]).*(I+1+[(I-1):-1:-I])),1);
%     Im=Ip';
%     Iph=kron(eye(2*J+1),Ip);
%     Imh=kron(eye(2*J+1),Im);
%     Ixh=(Iph+Imh)/2;
%     Iyh=(Iph-Imh)/2i;
    hJx=kron(Jx,Iid);
    hJy=kron(Jy,Iid);
    hJz=kron(Jz,Iid);
else
    hJx=Jx;
    hJy=Jy;
    hJz=Jz;
end

for ionn = 1:4
    [v,e]=eig(squeeze(Ham(ionn,:,:)));
    e=real(diag(e));
    e=e-min(e);
    [e,n]=sort(e);
    v=v(:,n);
    
    % calculate matric elements
    mJx = v'*hJx*v;
    mJy = v'*hJy*v;
    mJz = v'*hJz*v;
    % And corresponding weight factors:
    if t>0
        ne=exp(-e*11.6/t)/sum(exp(-e*11.6/t));
    else
        ne=(e==min(e)); %grundzustand (e(i)==min(e))=true=1 sonst 0
    end
%     erw=[sum(diag(mJx).*ne) sum(diag(mJy).*ne) sum(diag(mJz).*ne)];
%     mJx=mJx-erw(1)*eye(size(mJx));
%     mJy=mJy-erw(2)*eye(size(mJy));
%     mJz=mJz-erw(3)*eye(size(mJz));
    
    chi0 = double.empty(3,3,4,length(sample.omega),0);
    for n=1:length(sample.omega)
        w=(ones(size(v,1),1)*ne'-ne*ones(1,size(v,1)))./ ...
            (e*ones(1,size(v,1))-ones(size(v,1),1)*e'-sample.omega(n)-1i*sample.epsilon);
        
        chi0(:,:,ionn,n,1)=[sum(sum(mJx.'.*mJx.*w)) sum(sum(mJx.'.*mJy.*w)) sum(sum(mJx.'.*mJz.*w))
                            sum(sum(mJy.'.*mJx.*w)) sum(sum(mJy.'.*mJy.*w)) sum(sum(mJy.'.*mJz.*w))
                            sum(sum(mJz.'.*mJx.*w)) sum(sum(mJz.'.*mJy.*w)) sum(sum(mJz.'.*mJz.*w))];
    end
    chi0 = squeeze(chi0);
end

end
