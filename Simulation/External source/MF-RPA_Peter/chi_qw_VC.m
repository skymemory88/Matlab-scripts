function chi=chi_qw_VC(q,chi0,sample)
%berechnet chi(q,omega). q ist eine (anzahl_q mal 3) Matrix.

if strcmpi(sample.ionname,'Er')
    Jex = sample.exEr;          % Exchange coupling strength
    gL = sample.gLande_Er;      % Lande g-factor
    c = sample.ErHyp;           % fraction of ions with nuclear moment
end
if strcmpi(sample.ionname,'Ho')
    Jex = sample.exHo;          % Exchange coupling strength
    gL = sample.gLande_Ho;      % Lande g-factor
    c = sample.HoHyp;           % fraction of ions with nuclear moment
end


% Calculate dipole sum
N=4; % Number of magnetic atoms in unit cell
D=zeros(3,3,N,N,size(q,1));

for n=1:size(q,1)
    D(:,:,:,:,n)=dipole_direct(q(n,:),100)+exchange(q(n,:),Jex);
end

%%


% Convert from AA^-3 to meV (muB^2 =0.05368 meV AA^3)
%D=[-0.0301 0 0
%   0 -0.0301 0
%   0 0 0.0602];
vol=287.8917;
D=D*gL^2*0.05368;%*vol;


% % % % % % chi00 = chi0.hyp;
% % % % % % chi=zeros(3,3,size(chi00,4),size(q,1));
% % % % % % for nw=1:size(chi00,4)
% % % % % %     for nq=1:size(q,1)
% % % % % %         M=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
% % % % % %         for n=1:N  %n,m Summe ueber Ionen in der Einheitszelle
% % % % % %             for m=1:N
% % % % % %                 M((n-1)*3+(1:3),(m-1)*3+(1:3))=chi00(:,:,n,nw)*D(:,:,n,m,nq);
% % % % % %             end
% % % % % %         end
% % % % % %         chi(:,:,nw,nq)=1/4*[eye(3) eye(3) eye(3) eye(3)]*...
% % % % % %             ((eye(size(M))-M)\([chi00(:,:,1,nw);chi00(:,:,2,nw);chi00(:,:,3,nw);chi00(:,:,4,nw)]));
% % % % % %     end
% % % % % % end

if sample.hyp_nohyp
    chi0.hyp = chi0.nohyp;
    warning('Replace chi0(HF) with chi0(no HF)')
end
if sample.nohyp_hyp
    chi0.nohyp = chi0.hyp;
    warning('Replace chi0(no HF) with chi0(HF)')
end

%% Following the approach of J. Jensen, Rare-earth magnetism, p. 247
chi=zeros(3,3,length(sample.omega),size(q,1));

if ~sample.rmint
    for nw=1:length(sample.omega)
        for nq=1:size(q,1)
            M0=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
            M1=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
            
            for n=1:N  %n,m Summe ueber Ionen in der Einheitszelle
                for m=1:N
                    M1((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0.hyp(:,:,n,nw)*D(:,:,n,m,nq); %*(-1)^(n+m);
                    M0((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0.nohyp(:,:,n,nw)*D(:,:,n,m,nq); %*(-1)^(n+m);
                end
            end
            
            
            A  = eye(3*N) - c*M1 - (1-c)*M0;
            
            chi00 = 1/4*[eye(3) eye(3) eye(3) eye(3)]*(...
                            (1-c)*[ chi0.nohyp(:,:,1,nw);
                                    chi0.nohyp(:,:,2,nw);
                                    chi0.nohyp(:,:,3,nw);
                                    chi0.nohyp(:,:,4,nw)] ...
                + (M0/A)*(1-c)^2*[  chi0.nohyp(:,:,1,nw);
                                    chi0.nohyp(:,:,2,nw);
                                    chi0.nohyp(:,:,3,nw);
                                    chi0.nohyp(:,:,4,nw)]);
            
            chi10 = 1/4*[eye(3) eye(3) eye(3) eye(3)]*(...
                + (M1/A)*(1-c)*c*[  chi0.nohyp(:,:,1,nw);
                                    chi0.nohyp(:,:,2,nw);
                                    chi0.nohyp(:,:,3,nw);
                                    chi0.nohyp(:,:,4,nw)]);
            
            chi01 = 1/4*[eye(3) eye(3) eye(3) eye(3)]*(...
                + (M0/A)*(1-c)*c*[  chi0.hyp(:,:,1,nw);
                                    chi0.hyp(:,:,2,nw);
                                    chi0.hyp(:,:,3,nw);
                                    chi0.hyp(:,:,4,nw)]);
            
            chi11 = 1/4*[eye(3) eye(3) eye(3) eye(3)]*(...
                                c*[ chi0.hyp(:,:,1,nw);
                                    chi0.hyp(:,:,2,nw);
                                    chi0.hyp(:,:,3,nw);
                                    chi0.hyp(:,:,4,nw)] ...
                    + (M1/A)*c^2*[  chi0.hyp(:,:,1,nw);
                                    chi0.hyp(:,:,2,nw);
                                    chi0.hyp(:,:,3,nw);
                                    chi0.hyp(:,:,4,nw)]);
            
            
            
            chi(:,:,nw,nq) = chi00 + chi10 + chi01 + chi11;
            
            
        end
    end
else
    for nw=1:length(sample.omega)
        for nq=1:size(q,1)
            M0=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
            M1=zeros(3*N); %Blockmatrix, 3mal3 Bloecke entsprechen den D zwischen einzelnen Ionen in der Einehitszelle
            
            for n=1:N  %n,m Summe ueber Ionen in der Einheitszelle
                for m=1:N
                    M1((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0.hyp(:,:,n,nw)*D(:,:,n,m,nq); %*(-1)^(n+m);
                    M0((n-1)*3+(1:3),(m-1)*3+(1:3))=chi0.nohyp(:,:,n,nw)*D(:,:,n,m,nq); %*(-1)^(n+m);
                end
            end

            A  = eye(3*N) - c*M1 - (1-c)*M0;
            
            chi(:,:,nw,nq) = 1/4*[eye(3) eye(3) eye(3) eye(3)]/A*[eye(3); eye(3); eye(3); eye(3)];
        end
    end
    
end


return
