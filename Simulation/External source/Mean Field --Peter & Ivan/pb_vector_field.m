function pb_vector_field

clear all

opt = 2;

switch opt
    case 1
        sPHI = 90;
        Hz=0:1:8;

%         ll = 1;
%         for phi = sPHI
        [fields,altJ_H,J_H,nmax,ion] = LiYbF4_mf_sim(Hz,sPHI);
%             ll = ll + 1;
%         end

        % Make plots of Jx,Jy,Jz and staggered magnetisation
        option1(sPHI,Hz,fields,altJ_H,J_H,nmax,ion)
        
    case 2
        
%         sHx = [0   0.05 0.075 0.1 0.125 0.15 0.175 0.3];
%         sHy = 0.*sHx;
%         sHz = [0.3 0.3  0.3   0.3 0.3   0.3  0.3   0.3];
       

%         sHz = [0 0.1 0.2 0.3 0.4 0.5 0.6];
%         sHy = [0 0.1 0.2 0.3 0.4];
%         sHx = [0 0.05 0.1 0.15 0.2];

%         sHy = [0 0.05 0.1 0.15 0.2 0.3 0.4];
%         sHx = [0 0.05 0.1 0.15 0.2 0.3 0.4];
        sHz = [0 0.1 0.2 0.3 0.4 0.5 0.6];

        sHx = 0.*sHz;
        sHy = 0.*sHz;
       
%         sHx = [0.00];
%         sHy = [0.00];
%         sHz = [];

        for n = 1:length(sHx)
            Hx=sHx(n);
            Hy=sHy(n);
            Hz=sHz(n);

            ll = 1;
            % Create 3D model of the magnetic moments calculated by MF
            [fields{ll},altJ_H{ll},J_H{ll},nmax{ll},ion{ll}] = LiYbF4_mf_sim_2(Hx,Hy,Hz);

            % Look more carefully at what happens to the moments within the unit cell
            option2(fields,altJ_H,J_H,nmax,ion)
        end
end

end



function option1(sPHI,Hz,fields,altJ_H,J_H,nmax,ion)

%% Create plot ------------------------------------------------------------
nion=find(ion.prop);

hfig1 = figure(11);
clf
set(hfig1,'position',[10 10 600 600])
hold on
box on
h1=gca;

hfig2 = figure(12);
clf
set(hfig2,'position',[10 10 600 600])
hold on
box on
h2=gca;

hfig3 = figure(13);
clf
set(hfig3,'position',[10 10 600 600])
hold on
box on
h3=gca;

hfig4 = figure(14);
clf
set(hfig4,'position',[10 10 600 600])
hold on
box on
h4=gca;

hfig5 = figure(15);
clf
set(hfig5,'position',[10 10 600 600])
hold on
box on
h5=gca;

hfig6 = figure(16);
clf
set(hfig6,'position',[10 10 600 600])
hold on
box on
h6=gca;

% Hz=0:0.01:1;

ll = 1;
% for ll = 1:length(fields)
    p1(ll)=plot( h1,        Hz,altJ_H(1,:,nion));
    p2(ll)=plot( h2,        Hz,altJ_H(2,:,nion));
    p3(ll)=plot( h3,        Hz,altJ_H(3,:,nion));
    p4(ll)=plot( h4,        Hz,J_H(1,:,nion));
    p5(ll)=plot( h5,        Hz,J_H(2,:,nion));
    p6(ll)=plot( h6,        Hz,J_H(3,:,nion));
    set([p1(ll) p2(ll) p3(ll) p4(ll) p5(ll) p6(ll)],'color',setcolours((ll-1)/(length(fields)-1),'jet'))
% end

% for ll = 1:length(fields)
%     t.x= Hz;
%     t.y = altJ_H{ll}(2,:,nion);
%     t.e = 0.*Hz;
%     s(ll) = spec1d(t);
% end


%% Plot data
set([get(h1,'xlabel') get(h2,'xlabel') get(h3,'xlabel') get(h4,'xlabel') get(h5,'xlabel') get(h6,'xlabel')],...
    'string','H (T)')
ylabel(h1,'<(-1)^{i}J^x>')
ylabel(h2,'<(-1)^{i}J^y>')
ylabel(h3,'<(-1)^{i}J^z>')
ylabel(h4,'<J^x>')
ylabel(h5,'<J^y>')
ylabel(h6,'<J^z>')

% Save figures
saveplots(hfig1,'H-Jstag_x_0_-45to+45deg_v1')
saveplots(hfig2,'H-Jstag_y_0_-45to+45deg_v1')
saveplots(hfig3,'H-Jstag_z_0_-45to+45deg_v1')
saveplots(hfig4,'H-J_x_0_-45to+45deg_v1')
saveplots(hfig5,'H-J_y_0_-45to+45deg_v1')
saveplots(hfig6,'H-J_z_0_-45to+45deg_v1')

%% Analyse how Hc changes with applied field using staggered moment Jx
hfig7 = figure(17);
clf
set(hfig7,'position',[10 10 600 600])
hold on
box on
h7=gca;

for ll = 1:length(fields)
    t=altJ_H{ll}(1,:,nion)>1e-7;
    ind=find(t==1);
    
    Hc(ll) = Hz(ind(end));
end

p7(ll)=plot( h7,sPHI,Hc);
xlabel('\phi (deg)')
ylabel('H_c (T)')
saveplots(hfig7,'phi-Hc_-45to+45deg_v1')

%% Create colourmap of Jstag_y
hfig8 = figure(18);
clf
set(hfig8,'position',[10 10 600 600])
hold on
box on
h8=gca;

disp('I am pretty little girl')

mapplot(s,sPHI)
view(90,-90)
ylabel('H (T)')
xlabel('\phi (deg)')
colorbar
caxis([-1.5 1.5])
saveplots(hfig8,'Hc-phi_map_-45to+45deg_v1')

%% Convert the phase diagrams into Hx,0,Hz coordinates
Hx = [];
Hz = [];
sJx = [];
sJy = [];
sJz = [];
aJx = [];
aJy = [];
aJz = [];

for ll = 1:length(fields)
     Hx = [Hx; fields{ll}(1,:)'];
     Hz = [Hz; fields{ll}(3,:)'];
     sJx = [sJx; altJ_H{ll}(1,:,nion)'];
     sJy = [sJy; altJ_H{ll}(2,:,nion)'];
     sJz = [sJz; altJ_H{ll}(3,:,nion)'];
     
     aJx = [aJx; J_H{ll}(1,:,nion)'];
     aJy = [aJy; J_H{ll}(2,:,nion)'];
     aJz = [aJz; J_H{ll}(3,:,nion)'];     
end

[xx,yy] = meshgrid([-0.5:0.001:0.5],[0:0.001:1]);
[XX YY ZZ_sx] = griddata(Hx,Hz,sJx,xx,yy);
[XX YY ZZ_sy] = griddata(Hx,Hz,sJy,xx,yy);
[XX YY ZZ_sz] = griddata(Hx,Hz,sJz,xx,yy);
[XX YY ZZ_ax] = griddata(Hx,Hz,aJx,xx,yy);
[XX YY ZZ_ay] = griddata(Hx,Hz,aJy,xx,yy);
[XX YY ZZ_az] = griddata(Hx,Hz,aJz,xx,yy);


hfig21 = figure(21);
clf
set(hfig21,'position',[10 10 600 600])
hold on
box on
h21=gca;

hfig22 = figure(22);
clf
set(hfig22,'position',[10 10 600 600])
hold on
box on
h22=gca;

hfig23 = figure(23);
clf
set(hfig23,'position',[10 10 600 600])
hold on
box on
h23=gca;

hfig24 = figure(24);
clf
set(hfig24,'position',[10 10 600 600])
hold on
box on
h24=gca;

hfig25 = figure(25);
clf
set(hfig25,'position',[10 10 600 600])
hold on
box on
h25=gca;

hfig26 = figure(26);
clf
set(hfig26,'position',[10 10 600 600])
hold on
box on
h26=gca;

%
figure(hfig21)
p21=pcolor(h21,XX,YY,ZZ_sx);
hc=colorbar;
ylabel(hc,'<(-1)^iJ^x>')
grid on
set(gca,'xtick',[-0.5:0.1:0.5])
%
figure(hfig22)
p22=pcolor(h22,XX,YY,ZZ_sy);
hc=colorbar;
ylabel(hc,'<(-1)^iJ^y>')
grid on
caxis([-1.6 1.6])
set(gca,'xtick',[-0.5:0.1:0.5])
%
figure(hfig23)
p23=pcolor(h23,XX,YY,ZZ_sz);
hc=colorbar;
ylabel(hc,'<(-1)^iJ^z>')
grid on
set(gca,'xtick',[-0.5:0.1:0.5])
%
figure(hfig24)
p24=pcolor(h24,XX,YY,ZZ_ax);
hc=colorbar;
ylabel(hc,'<J^x>')
grid on
caxis([-1.6 1.6])
set(gca,'xtick',[-0.5:0.1:0.5])
%
figure(hfig25)
p25=pcolor(h25,XX,YY,ZZ_ay);
hc=colorbar;
ylabel(hc,'<J^y>')
grid on
set(gca,'xtick',[-0.5:0.1:0.5])
%
figure(hfig26)
p26=pcolor(h26,XX,YY,ZZ_az);
hc=colorbar;
ylabel(hc,'<J^z>')
grid on
set(gca,'xtick',[-0.5:0.1:0.5])

set([p21 p22 p23 p24 p25 p26],'edgecolor','none')
set([get(h21,'xlabel') get(h22,'xlabel') get(h23,'xlabel') get(h24,'xlabel') get(h25,'xlabel') get(h26,'xlabel')],...
    'string','H_x (T)')
set([get(h21,'ylabel') get(h22,'ylabel') get(h23,'ylabel') get(h24,'ylabel') get(h25,'ylabel') get(h26,'ylabel')],...
    'string','H_z (T)')

saveplots(hfig21,'H-Jstag_x_map_-45to+45deg_v1')
saveplots(hfig22,'H-Jstag_y_map_-45to+45deg_v1')
saveplots(hfig23,'H-Jstag_z_map_-45to+45deg_v1')
saveplots(hfig24,'H-J_x_map_-45to+45deg_v1')
saveplots(hfig25,'H-J_y_map_-45to+45deg_v1')
saveplots(hfig26,'H-J_z_map_-45to+45deg_v1')

end

function option2(fields,altJ_H,J_H,nmax,ion)

nion=find(ion{end}.prop);

curdir = cd;
cd('C:\Users\ikovacev\Desktop\LiYbF4_MF-calc')

lattice.a = ion{1}.a{nion}(1,1);
lattice.b = ion{1}.a{nion}(2,2);
lattice.c = ion{1}.a{nion}(3,3);
lattice = basis(lattice);

r_pos = [   [1/2 0 0]
            [1/2 1/2 1/4]
            [0 1/2 1/2]
            [0 0 3/4]    ]*ion{1}.a{nion};


        
        
for nfields = 1:length(fields)

    
    hfig = figure(200);
    clf
    hold on
    axis equal
    axis off
    xlim([-1 6])
    ylim([-1 6])
    zlim([-2 12])
    
    l = drawunitcell(lattice);
    draw_coord_system(lattice)
    
    set(gca,'DataAspectRatio',[1 1 1],...
            'PlotBoxAspectRatio',[1 1 1]);
    for natoms = 1:4
        disp([num2str(natoms) '. ' num2str(norm(ion{nfields}.mom(natoms,:,nion)),' %3.3f')])
        [hline(:,natoms),hhead(:,natoms)] = arrow3d(r_pos(natoms,:),r_pos(natoms,:)+ion{nfields}.mom(natoms,:,nion));
    end
    set(hline,'facecolor',[0.2 0.2 0.7])
    set(hhead,'facecolor',[0.2 0.2 0.7])
    camlight headlight
    
    m1 = ['m_1 = [' num2str(ion{nfields}.mom(1,:,nion),' %3.2f') '], (' num2str(norm(ion{nfields}.mom(1,:,nion)),'%3.2f') ')'];
    m2 = ['m_2 = [' num2str(ion{nfields}.mom(2,:,nion),' %3.2f') '], (' num2str(norm(ion{nfields}.mom(2,:,nion)),'%3.2f') ')'];
    m3 = ['m_3 = [' num2str(ion{nfields}.mom(3,:,nion),' %3.2f') '], (' num2str(norm(ion{nfields}.mom(3,:,nion)),'%3.2f') ')'];
    m4 = ['m_4 = [' num2str(ion{nfields}.mom(4,:,nion),' %3.2f') '], (' num2str(norm(ion{nfields}.mom(4,:,nion)),'%3.2f') ')'];
    
    % Calculate the Zeeman interaction between moments and external
    % magnetic field
    DE = -sum(dot(ion{nfields}.mom(:,:,nion),repmat(fields{nfields}',4,1),2));
    
    title(['H = [' num2str(fields{nfields}',' %3.3f') '], U_Z = ' num2str(DE, ' %3.3f')])
    t = annotation('textbox',[0.01 0.1 0.4 0.1],'string',char(m1,m2,m3,m4),'edgeColor','none');
    
    % Orthographic view
    view(130,50)
%     pause(1)
    saveplots(hfig,['H = [' num2str(fields{nfields}',' %3.3f') '] 3D'])
    save2latex(hfig,['H_' num2str(fields{nfields}(1)*1000,' %d') '_' num2str(fields{nfields}(2)*1000,' %d') '_' num2str(fields{nfields}(3)*1000,' %d') '_3D'])
    delete(t)
    title('')
    
%     % View along the -c-axis
%     view(90,90)
% %     pause(1)
%     saveplots(hfig,['H = [' num2str(fields{nfields}',' %3.3f') '] ab-plane'])
%     save2latex(hfig,['H_' num2str(fields{nfields}(1)*1000,' %d') '_' num2str(fields{nfields}(3)*1000,' %d') '_ab-plane'])
%     
%     
%     % View along the -a-axis
%     view(90,0)
% %     pause(1)
%     saveplots(hfig,['H = [' num2str(fields{nfields}',' %3.3f') '] bc-plane'])
%     save2latex(hfig,['H_' num2str(fields{nfields}(1)*1000,' %d') '_' num2str(fields{nfields}(3)*1000,' %d') '_bc-plane'])
%     
%     % View along the -b-axis
%     view(-180,0)
% %     pause(1)
%     saveplots(hfig,['H = [' num2str(fields{nfields}',' %3.3f') '] ac-plane'])
%     save2latex(hfig,['H_' num2str(fields{nfields}(1)*1000,' %d') '_' num2str(fields{nfields}(3)*1000,' %d') '_ac-plane'])
    
end


    
cd(curdir)

end

%% Sub-functions ----------------------------------------------------------

function [fields,altJ_H,J_H,nmax,ion] = LiYbF4_mf_sim(Hz,phi)


%% Setup parameters.
global strategies; % global convergence strategies switches
strategies.powerlaw=false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=false; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.


global rundipole
rundipole = true;

%<J> as a function of H (Q=1) or T (Q=0) ?
Q=1;


if Q
    phi = phi*pi/180;                                 % rotation angle          
    
    Rx = [1         0           0
          0         cos(phi)    sin(phi);
          0        -sin(phi)    cos(phi)];          % rotation about z 
      
    Ry = [cos(phi)  0          -sin(phi);
          0         1           0;
          sin(phi)  0           cos(phi)];          % rotation about z 
      
    Rz = [cos(phi)  sin(phi)    0;
         -sin(phi)  cos(phi)    0;
          0         0           1];                 % rotation about z

     
    temp=0;
%     Hz=[0:0.4:1.2 1.25:0.05:1.6];
%     Hz=0:0.01:1;
%     Hz=2:-0.005:0;
%     Hz = [0:1:15 15.1:0.1:20 20.01:0.01:25 25:0.1:30];
    fields=Ry*[zeros(size(Hz))
            zeros(size(Hz))
            Hz
                       ];
else
    temp=0:0.001:0.5; %  range of temperatures
    fields=[0;0;0];
end



%Ions' names
ion.name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
ion.name_hyp=[{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];

%Ions' proportions
ion.prop=[0;1;0;0;0;0];

%Ions' J, L, S values
%      Er       Ho      Yb      Tm      Gd      Y
ion.J=[15/2;    8;      7/2;    6;      7/2;    1];
ion.L=[6;       6;      3;      5;      0;      1];
ion.S=[3/2;     2;      1/2;    1;      7/2;    1];

%Ions' lattice parameters
ion.a=[{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
       {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
       {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
       {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y

%Ions' cf parameters
% ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
%        [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
%        [663   12.5   102   -0.62   -16   0]                             %Yb        
%        [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
%        [0 0 0 0 0 0]                                                    %Gd
%        [0 0 0 0 0 0]]/1000;                                             %Y

% Ions' cf parameters
ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
       [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
       [646.2016   15.3409  116.4854   -0.0686  -15.1817    0.0000]     %Yb        
       [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
       [0 0 0 0 0 0]                                                    %Gd
       [0 0 0 0 0 0]]/1000;                                             %Y
 
   
%Ions' renorm
ion.renorm=[[1,1,1]         %Er
            [1,1,0.785]     %Ho
            [1,1,1]         %Yb
            [1,1,1]         %Tm
            [1,1,1]         %Gd
            [1,1,1]];       %Y

%Ions' hyperfine
ion.hyp=[0;1;0;0;0;0];

%Ions' nuclear spin I
ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[0.00043412; 0.003361; 0; 0; 0; 0];

%Exchange
ion.ex=[0; -0.0001; 0; 0; 0; 0];
% ex.Er=0; % magnetic exchange in J*S_iS_j
% ex.Ho=-0.0001; %0.1 microeV from Henrik and Jens PRB
% ex.Yb=0;
% ex.Tm=0;

% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

%Ions' initial moments
ion.mom(:,:,1)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Er
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

demagn=true;

alpha=1; % shape is a sphere.

%% Calculate

% profile on
% profile clear
[ion,history,E]=LiIonsF4(ion,temp,fields,demagn,alpha);
% profile report

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result

for m = 1:length(history)
    niter(m) = history{m}.iterations;
end
nmax = max(niter);
    
if Q
    altJ_H=zeros([3,size(fields,2),size(ion.name,1)]);
    J_H=zeros([3,size(fields,2),size(ion.name,1)]);  
    Jnorm_H=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_H(:,:,i)=squeeze(mean(ion.altJs(:,:,:,:,i),1));
        J_H(:,:,i)=squeeze(mean(ion.Js(:,:,:,:,i),1));
        Jnorm_H(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
else
    altJ_T=zeros([3,size(temp,2),size(ion.name,1)]);
    J_T=zeros([3,size(temp,2),size(ion.name,1)]); 
    Jnorm_T=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_T(:,:,i)=mean(ion.altJs(:,:,:,i),1);
        J_T(:,:,i)=mean(ion.Js(:,:,:,i),1);
        Jnorm_T(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
end

% nion=find(ion.prop);


% hfig1 = figure(1);
% clf
% set(hfig1,'position',[10 10 600 600])
% 
% if Q
% %     p=plot(Hz,Jnorm_H(:,1));
% %     p=plot(Hz,altJ_H(1,:,n),Hz,J_H(3,:,n));
% % %     legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>');
% %     legend('<(-1)^{i}J^x>','<J^z>');
% %     title('LiYbF_4')
% % %     title(sprintf('LiGdF_4'))
% %     xlabel('H_z [T]');
% %     ylabel('|J|');
% %     saveplots(hfig1,['Simulation of LiYbF4 Hz vs J'])
% 
%     p=plot(Hz,altJ_H(1,:,n),Hz,altJ_H(2,:,n),Hz,altJ_H(3,:,n),Hz,J_H(1,:,n),Hz,J_H(2,:,n),Hz,J_H(3,:,n));
% %     legend('<(-1)^{i}J^x_{Er}>','<J^z_{Er}>');
%     legend('<(-1)^{i}J^x>','<(-1)^{i}J^y>','<(-1)^{i}J^z>','<J^x>','<J^y>','<J^z>');
%     title(char('LiYbF_4',...
%             ['B = ' num2str(ion.B(3,:)*1000) ' meV'],...
%             ['H = H_0(' num2str(-sin(phi),3) ', 0, ' num2str(cos(phi),3) '), \phi(deg) = ' num2str(phi*180/pi)]))
%     
% %     title(sprintf('LiGdF_4'))
%     xlabel('H_z [T]');
% else
%     p=plot(temp,altJ_T(1,:,n),temp,J_T(3,:,n));
%     legend('<(-1)^{i}J^x>','<J^z>');
% %     title(sprintf('LiEr_{%1.1f}Ho_{%1.1f}Yb_{%1.1f}Tm_{%1.1f}Gd_{%1.1f}Y_{%1.1f}F_4',ion.prop(1),ion.prop(2),ion.prop(3),ion.prop(4),ion.prop(5),ion.prop(6)))
%     title('LiYbF_4')
%     xlabel('T [K]')
% %     ylabel('|J|');
% %     saveplots(hfig1,['Simulation of LiYbF4 T vs J'])
% end

% print -dpng example.png

%% Plot history
% cols=jet(length(temp)-1);
% iterations=0;
% figure(2)
% for t=2:length(temp)
% %    relax=squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2);
%     relax=sqrt((squeeze(mean(history{t}.ho(:,:,3),2))-mean(history{t}.ho(end,:,3),2)).^2+(squeeze(history{t}.er(:,1,1))-history{t}.er(end,1,1)).^2);
%     line(temp(t)*ones(size(relax)),1:size(relax,1),relax,'Color',cols(t-1,:));
%     hold on
%     iterations=iterations+size(relax,1);
% end
% fprintf('iterations=%d\n',iterations);
% set(gca,'Box','on')
% set(gca,'XGrid','on')
% set(gca,'YGrid','on')
% set(gca,'ZGrid','on')
% xlabel('Temperature');
% ylabel('Iterations');
% zlabel('Ho J_i^z-<J^z>');


end

function [fields,altJ_H,J_H,nmax,ion] = LiYbF4_mf_sim_2(Hx,Hy,Hz)


%% Setup parameters.
global strategies; % global convergence strategies switches
strategies.powerlaw=false; % Use a starting guess with the assumption the variables follow a power law.
strategies.accelerator=0.0; % accelerate. Beware: it often doesn't accelerate actually. Can speed up if the system has to break a symetry
strategies.damping=0.1; % damping factor. May reduce exponential fit efficiency
strategies.expfit=false; % turn on fitting to an exponential.
strategies.expfit_period=50; % period between exponential fit.
strategies.expfit_deltaN=20; % The exponential fit points distance. Three points used: N, N-deltaN, N-2deltaN.


global rundipole
rundipole = true;

%<J> as a function of H (Q=1) or T (Q=0) ?
Q=1;


if Q
%     phi = phi*pi/180;                                 % rotation angle          
%     
%     Rx = [1         0           0
%           0         cos(phi)    sin(phi);
%           0        -sin(phi)    cos(phi)];          % rotation about z 
%       
%     Ry = [cos(phi)  0          -sin(phi);
%           0         1           0;
%           sin(phi)  0           cos(phi)];          % rotation about z 
%       
%     Rz = [cos(phi)  sin(phi)    0;
%          -sin(phi)  cos(phi)    0;
%           0         0           1];                 % rotation about z

     
    temp=0;
%     Hz=[0:0.4:1.2 1.25:0.05:1.6];
%     Hz=0:0.01:1;
%     Hz=2:-0.005:0;
%     Hz = [0:1:15 15.1:0.1:20 20.01:0.01:25 25:0.1:30];
%     fields=Ry*[zeros(size(Hz))
%             zeros(size(Hz))
%             Hz
%                        ];

    fields = [Hx; Hy; Hz];
else
    temp=0:0.001:0.5; %  range of temperatures
    fields=[0;0;0];
end



%Ions' names
ion.name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];
ion.name_hyp=[{'Er_hyp'};{'Ho_hyp'};{'Yb_hyp'};{'Tm_hyp'};{'Gd_hyp'};{'Y_hyp'}];

%Ions' proportions
ion.prop=[1;0;0;0;0;0];

%Ions' J, L, S values
%      Er       Ho      Yb      Tm      Gd      Y
ion.J=[15/2;    8;      7/2;    6;      7/2;    1];
ion.L=[6;       6;      3;      5;      0;      1];
ion.S=[3/2;     2;      1/2;    1;      7/2;    1];

%Ions' lattice parameters
ion.a=[{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
       {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
       {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
       {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gdk
       {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y

%Ions' cf parameters
% ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
%        [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
%        [663   12.5   102   -0.62   -16   0]                             %Yb        
%        [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
%        [0 0 0 0 0 0]                                                    %Gd
%        [0 0 0 0 0 0]]/1000;                                             %Y

% Ions' cf parameters
ion.B=[[60.2258   -0.1164   -4.3280   -0.0019   -0.0850   -0.0227]      %Er
       [-60   0.35   3.6   0.0004   0.07   0.006]                       %Ho
       [646.2016   15.3409  116.4854   -0.0686  -15.1817    0.0000]     %Yb        
       [224.3   -1.85   -11.7   0.002   0.2645   0.1377]                %Tm
       [0 0 0 0 0 0]                                                    %Gd
       [0 0 0 0 0 0]]/1000;                                             %Y
 
   
%Ions' renorm
ion.renorm=[[1,1,1]         %Er
            [1,1,0.785]     %Ho
            [1,1,1]         %Yb
            [1,1,1]         %Tm
            [1,1,1]         %Gd
            [1,1,1]];       %Y

%Ions' hyperfine
ion.hyp=[0;0;0;0;0;0];

%Ions' nuclear spin I
ion.I=[3; 3.5; 0; 0; 0; 0];
ion.A=[0.00043412; 0.003361; 0; 0; 0; 0];

%Exchange
ion.ex=[0; -0.0001; 0; 0; 0; 0];
% ex.Er=0; % magnetic exchange in J*S_iS_j
% ex.Ho=-0.0001; %0.1 microeV from Henrik and Jens PRB
% ex.Yb=0;
% ex.Tm=0;

% exHo=-0.000542; % from Bitko et al.
% exHo=-0.000436; % from Conradin ?

%Ions' initial moments
% ion.mom(:,:,1)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Er BLAFM x !!!
% ion.mom(:,:,1)=[0 1 0; 0 1 0; 0 -1 0; 0 -1 0];      %Er BLAFM y !!!
ion.mom(:,:,1)=[-1 1 0; 1 1 0; 1 -1 0; -1 -1 0];      %Er spiral xy 45d clockwise  !!!
% ion.mom(:,:,1)=[-1 -1 0; 1 -1 0; 1 1 0; -1 1 0];      %Er spiral xy 45d anticlockwise !!! 
% ion.mom(:,:,1)=[1 1 0; -1 1 0; -1 -1 0; 1 -1 0];      %Er spiral xy 45d any linear combination of x BLAFM and y BLAFM  !!!
ion.mom(:,:,2)=[0 0 1;  0 0 1;  0 0 1; 0 0 1]*2.6;  %Ho 
ion.mom(:,:,3)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Yb BLAFM x !!!
% ion.mom(:,:,3)=[0 -1 0; 0 -1 0; 0 1 0; 0 1 0];      %Yb BLAFM y !!!
% ion.mom(:,:,3)=[1 0 0;  0 1 0; -1 0 0; 0 -1 0];      %Yb spiral clockwise
% -> +45 degrees rotation 
% ion.mom(:,:,3)=[1 0 0;  0 -1 0; -1 0 0; 0 1 0];      %Yb spiral
% anti-clockwise -> -45 degrees rotation 
% ion.mom(:,:,3)=[1 0 0;  -1 0 0; 0 1 0; 0 -1 0];      %Yb xy AFM -> BLAFM x
% ion.mom(:,:,3)=[0 1 0; -1 0 0; 1 0 0;  0 -1 0];      %Yb xy AFM -> BLAFM y
% ion.mom(:,:,3)=[1 0 0; 1 0 0; -1 0 0; -1 0 0];      %Yb BLAFM x -> BLAFM y
% ion.mom(:,:,3)=[0 -1 0; 0 1 0; 0 1 0; 0 -1 0];      %Yb BLAFM y -> BLAFM x
% ion.mom(:,:,3)=[-1 1 0; 1 1 0; 1 -1 0; -1 -1 0];      %Yb spiral xy 45d anticlockwise  !!!
% ion.mom(:,:,3)=[1 1 0; -1 1 0; -1 -1 0; 1 -1 0];      %Yb spiral xy 45d any linear combination of x BLAFM and y BLAFM  !!!
% ion.mom(:,:,3)=[-1 -1 0; 1 -1 0; 1 1 0; -1 1 0];      %Yb spiral xy 45d clockwise !!! 
% ion.mom(:,:,3)=[-1 -1 0; -1 -1 0; 1 1 0; 1 1 0];      %Yb BLAFM 45d -> BLAFM y
% ion.mom(:,:,3)=[1 1 0; -1 -1 0; -1 -1 0; 1 1 0];      %Yb BLAFM 45d -> BLAFM x
% ion.mom(:,:,3)=[1 1 0; -1 -1 0; -1 1 0; -1 1 0];      %Yb xy AFM 45d -> anticlockwise 45d
ion.mom(:,:,4)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Tm
ion.mom(:,:,5)=[1 0 0; -1 0 0; -1 0 0; 1 0 0];      %Gd
ion.mom(:,:,6)=[0 0 0; 0 0 0; 0 0 0; 0 0 0];      %Y

demagn=true;

alpha=1; % shape is a sphere.

%% Calculate

% profile on
% profile clear
[ion,history,E]=LiIonsF4(ion,temp,fields,demagn,alpha);
% profile report

%JEr is the mean <J>. altJEr is the alternated mean <(-1)^(i+j)*J_(i,j)>

    %% Plot result

for m = 1:length(history)
    niter(m) = history{m}.iterations;
end
nmax = max(niter);
    
if Q
    altJ_H=zeros([3,size(fields,2),size(ion.name,1)]);
    J_H=zeros([3,size(fields,2),size(ion.name,1)]);  
    Jnorm_H=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_H(:,:,i)=squeeze(mean(ion.altJs(:,:,:,:,i),1));
        J_H(:,:,i)=squeeze(mean(ion.Js(:,:,:,:,i),1));
        Jnorm_H(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
else
    altJ_T=zeros([3,size(temp,2),size(ion.name,1)]);
    J_T=zeros([3,size(temp,2),size(ion.name,1)]); 
    Jnorm_T=zeros([size(fields,2),size(ion.name,1)]);
    for i=1:size(ion.name,1)
        altJ_T(:,:,i)=mean(ion.altJs(:,:,:,i),1);
        J_T(:,:,i)=mean(ion.Js(:,:,:,i),1);
        Jnorm_T(:,i)=squeeze(mean(ion.Jmom_norm(:,:,i),1));
    end
end




end

function col = setcolours(perc,type)
% Function to generate a colour scheme for plots
%
% INPUT:
% perc                  fraction of the colour vector to use
% type                  different types of plot colour schemes
%
% OUTPUT:
% col                   1x3 vector containing rgb colour information

try
    cmap = colormap(type);
catch ME
    type = 'jet';
    cmap = colormap(type);
end


switch lower(type)
        
    case 'rand'
        col = rand(1,3);
        
    otherwise
        % Load the specified colour map
        N = length(cmap);
        col = interp1(cmap,perc*(N-1)+1 );

end

end

function saveplots(hfig,figname)

curdir = cd;
cd('C:\Users\ikovacev\Desktop\LiYbF4_MF-calc\Figures')

% saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
% print2eps(figname,hfig)
% [~,~] = eps2xxx([figname '.eps'],{'jpeg','pdf'});

disp(['Figure ' figname ' saved to '])
disp(cd)
cd(curdir)

end

function save2latex(hfig,figname)

curdir = cd;
cd('C:\Users\ikovacev\Desktop\LiYbF4_MF-calc\Figures')

print2eps(figname,hfig)
cd(curdir)

end

function l = drawunitcell(lattice)
% Draw outline of the unit cell

    options.unitcell.cyl_vs_lin = 'lin';
    options.unitcell.lin_style = '-';
    options.unitcell.lin_wid = 1;
    options.unitcell.col = 'k';
    
    U{1} = [0 0 0];
    U{2} = [1 0 0];
    U{3} = [1 1 0];
    U{4} = [0 1 0];
    
    U{5} = [0 0 1];
    U{6} = [1 0 1];
    U{7} = [1 1 1];
    U{8} = [0 1 1];
    
    % for m = 1:8
    %     U{m} = a*U{m}(1) + b*U{m}(2) + c*U{m}(3) + r{1};
    % end
    
    for m = 1:8
        U{m} = lattice.Ra*U{m}(1) + lattice.Rb*U{m}(2) + lattice.Rc*U{m}(3);
    end
    
    
    
    switch lower(options.unitcell.cyl_vs_lin)
        
        case 'cyl'
            
            options.unitcell.lin_wid = 0.001*options.unitcell.lin_wid;

            drawcylinder(U{1},U{2},options)
            drawcylinder(U{2},U{3},options)
            drawcylinder(U{3},U{4},options)
            drawcylinder(U{4},U{1},options)

            drawcylinder(U{1},U{5},options)
            drawcylinder(U{2},U{6},options)
            drawcylinder(U{3},U{7},options)
            drawcylinder(U{4},U{8},options)

            drawcylinder(U{5},U{6},options)
            drawcylinder(U{6},U{7},options)
            drawcylinder(U{7},U{8},options)
            drawcylinder(U{8},U{5},options)

        case 'lin'
            
            l(1) =  line([U{1}(1) U{2}(1)],[U{1}(2) U{2}(2)],[U{1}(3) U{2}(3)]);
            l(2) =  line([U{2}(1) U{3}(1)],[U{2}(2) U{3}(2)],[U{2}(3) U{3}(3)]);
            l(3) =  line([U{3}(1) U{4}(1)],[U{3}(2) U{4}(2)],[U{3}(3) U{4}(3)]);
            l(4) =  line([U{4}(1) U{1}(1)],[U{4}(2) U{1}(2)],[U{4}(3) U{1}(3)]);

            l(5) =  line([U{1}(1) U{5}(1)],[U{1}(2) U{5}(2)],[U{1}(3) U{5}(3)]);
            l(6) =  line([U{2}(1) U{6}(1)],[U{2}(2) U{6}(2)],[U{2}(3) U{6}(3)]);
            l(7) =  line([U{3}(1) U{7}(1)],[U{3}(2) U{7}(2)],[U{3}(3) U{7}(3)]);
            l(8) =  line([U{4}(1) U{8}(1)],[U{4}(2) U{8}(2)],[U{4}(3) U{8}(3)]);

            l(9) =  line([U{5}(1) U{6}(1)],[U{5}(2) U{6}(2)],[U{5}(3) U{6}(3)]);
            l(10)=  line([U{6}(1) U{7}(1)],[U{6}(2) U{7}(2)],[U{6}(3) U{7}(3)]);
            l(11)=  line([U{7}(1) U{8}(1)],[U{7}(2) U{8}(2)],[U{7}(3) U{8}(3)]);
            l(12)=  line([U{8}(1) U{5}(1)],[U{8}(2) U{5}(2)],[U{8}(3) U{5}(3)]);

            set(l,...
                'Marker','none',...
                'LineStyle',options.unitcell.lin_style,...
                'linewidth',options.unitcell.lin_wid,...
                'color',options.unitcell.col)

    end

end

function lattice = basis(lattice)
% Calculate the basis vectors of the unit cell, {Ra, Rb, Rc}

if ~isfield(lattice,'al')
    lattice.al = 90;
end
if ~isfield(lattice,'be')
    lattice.be = 90;
end
if ~isfield(lattice,'ga')
    lattice.ga = 90;
end

al = lattice.al*pi/180;
be = lattice.be*pi/180;
ga = lattice.ga*pi/180;

lattice.Ra = [lattice.a 0 0];
lattice.Rb = lattice.b*[cos(ga) sin(ga) 0];

lattice.Rc = [0 0 0];
lattice.Rc(1) = lattice.c*cos(be);
lattice.Rc(2) = lattice.c*(cos(al) - cos(be)*cos(ga))/sin(ga);
if norm(lattice.Rc(2)) > lattice.Rc; error('Angles are geometrically inconsistent'), end
lattice.Rc(3) = sqrt(lattice.c^2 - lattice.Rc(1)^2 - lattice.Rc(2)^2);


end

function draw_coord_system(lattice)
%% DRAW BASIS VECTORS OF THE LATTICE:
% Draw arrows representing the coordinate system
options.coords.plot = 'n';              % Draw the coordinate system y/[n]
options.coords.r = -[1 1 1]*0.5;             % Position of the axes
options.coords.len = 1;                 % Length of the arrows
options.coords.norm_vs_act = 'norm';    % Normlise vector lengths
options.coords.col = [0 0 0];           % Colour of the axis arrows
options.coords.ax_x_labels = '$a$';     % Label the x,y,z axis (use x,y,z
                                        % by default)
options.coords.ax_y_labels = '$b$';       
options.coords.ax_z_labels = '$c$';
%                                           Position of the axes labels:
options.coords.ax_x_r = options.coords.r + [0.1 0 -0.1]*options.coords.len;      
options.coords.ax_y_r = options.coords.r + [0.1 0.1 0]*options.coords.len;
options.coords.ax_z_r = options.coords.r + [0.1 0 0.1]*options.coords.len;
options.coords.ax_ftsz = 20;            % Font size of the labels
options.coords.arrow_ang = 20;          % Definitions for the arrows
options.coords.arrow_p = [0.25,0.2];    % p(1) is the ratio of arrow head height 
                                        % to the distance between start and 
                                        % stop points P(2) is the ratio of 
                                        % arrow body
options.coords.arrow_n = [20,10];       % The arrow head has n(1) equally 
                                        % spaced points around its circumference.
                                        % The arrow body has n(2) equally spaced 
                                        % points around its circumference.
% -------------------------------------------------------------------------



% draw arrows to represent the axes

col = options.coords.col;
p = options.coords.arrow_p;
ang = options.coords.arrow_ang;
n = options.coords.arrow_n;


R0 = options.coords.r;                      % origin of coordinate axes

ra = lattice.Ra;
rb = lattice.Rb;
rc = lattice.Rc;

if strcmpi(options.coords.norm_vs_act,'norm')
    ra = ra/norm(ra);
    rb = rb/norm(rb);
    rc = rc/norm(rc);
end

ra = ra*options.coords.len;
rb = rb*options.coords.len;
rc = rc*options.coords.len;

% x-axis
[h1(1:2) h2(1:2)] = arrow3d(R0,R0 + ra,ang,'cylinder',p,n);
text(   R0(1) + ra(1) +  options.coords.ax_x_r(1),...
        R0(2) + ra(2) +  options.coords.ax_x_r(2),...
        R0(3) + ra(3) +  options.coords.ax_x_r(3),...
                options.coords.ax_x_labels,...
                'fontsize',options.coords.ax_ftsz,...
                'interpreter','latex')

% y-axis
[h1(3:4) h2(3:4)] = arrow3d(R0,R0 + rb,ang,'cylinder',p,n);
text(   R0(1) + rb(1) +  options.coords.ax_y_r(1),...
        R0(2) + rb(2) +  options.coords.ax_y_r(2),...
        R0(3) + rb(3) +  options.coords.ax_y_r(3),...
                options.coords.ax_y_labels,...
                'fontsize',options.coords.ax_ftsz,...
                'interpreter','latex')

% z-axis
[h1(5:6) h2(5:6)] = arrow3d(R0,R0 + rc,ang,'cylinder',p,n);
text(   R0(1) + rc(1) +  options.coords.ax_z_r(1),...
        R0(2) + rc(2) +  options.coords.ax_z_r(2),...
        R0(3) + rc(3) +  options.coords.ax_z_r(3),...
                options.coords.ax_z_labels,...
                'fontsize',options.coords.ax_ftsz,...
                'interpreter','latex')


set(h1,'FaceColor',col)
set(h2,'FaceColor',col)

end
