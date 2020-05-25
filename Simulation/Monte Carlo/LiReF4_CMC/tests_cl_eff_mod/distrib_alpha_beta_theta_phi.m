N=50000;
moms=zeros(N,3);
theta=zeros(N,1);
phi=zeros(N,1);
alpha=zeros(N,1);
beta=zeros(N,1);

%Ions' names
name=[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}];

%Ions' J, L, S values
J=[15/2; 8; 7/2; 6; 7/2; 1];
L=[6; 6; 3; 5; 0; 1];
S=[3/2; 2; 1/2; 1; 7/2; 1];

%Ions' lattice parameters
a=[{[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Er
   {[5.175 0 0; 0 5.175 0; 0 0 10.75]}      %Ho
   {[5.132 0 0; 0 5.132 0; 0 0 10.59]}      %Yb
   {[5.150 0 0; 0 5.150 0; 0 0 10.64]}      %Tm
   {[5.162 0 0; 0 5.162 0; 0 0 10.70]}      %Gd
   {[5.132 0 0; 0 5.132 0; 0 0 10.59]}];    %Y

%Ions' cf parameters
B=[[60.2258   -0.1164   -4.3280   0  -0.0019   -0.0850   -0.0227]     %Er [60.2258   -0.1164   -4.3280   0  -0.0019   -0.0850   -0.0227]  [63   -0.55   -5.54  0.47  -0.000006   -0.1082   -0.0146]
   [-60   0.35   3.6  0  0.0004   0.07   0.006]                       %Ho
   [663   12.5   102  0  -0.62   -16   0]                             %Yb !!! Conradin's parameters !!!
   [224.3   -1.85   -11.7  -15.2  0.002   0.2645   0.1377]            %Tm
   [0 0 0 0 0 0 0]                                                    %Gd
   [0 0 0 0 0 0 0]]/1000;                                             %Y
    
%B(1,:)=[60.2258 0 0 0 0 0 0]/1000;

F=[0.005; 0; 0; 0; 0; 0]; %cf Science LiErF4 2012, supp mat
params.field_changed=1;

for i=3
    ion{i}.name=name{i};
    ion{i}.num=i;
    ion{i}.J=J(i);
    ion{i}.L=L(i);
    ion{i}.S=S(i);
    ion{i}.B=B(i,:);  
    ion{i}.F=F(i);   
    ion{i}.gLande=gLande(ion{i}.L,ion{i}.S); 
    ion{i}.Hcf=Hnew;
    [ion{i}.VV,~]=Ising_basis([0,0,0],ion{i});
end 
        
%%

colnames='angles0pi angles02pi alpha beta theta phi';
angles0pi=0:pi/200:pi;
angles02pi=0:pi/100:2*pi;
n=length(angles0pi);

for j=3
    for t=1:N
        [alpha(t),beta(t)]=random_angles;
        [moms(t,:),~,params]=cl_eff_mod(alpha(t),beta(t),[0,0,0],ion{j},params);
        [phi(t),theta(t),~]=cart2sph(moms(t,1),moms(t,2),moms(t,3));
    end
end
    


    %%
    figure(70*j+1)
    hist(alpha,500)
    xlabel('\alpha')
    titre=sprintf('For %d rotations, Li%sF_4',N,ion{j}.name);
    title(titre)

    figure(90*j+2)
    hist(beta,500)
    xlabel('\beta')
    titre=sprintf('For %d rotations, Li%sF_4',N,ion{j}.name);
    title(titre)
    
    figure(500*j+3)
    hist(theta,500)
    xlabel('\theta')
    titre=sprintf('Distribution of the directions of the moments for %d rotations, Li%sF_4',N,ion{j}.name);
    title(titre)
   
    figure(10*j+4)
    hist(phi,500)
    xlabel('\phi')
    titre=sprintf('Distribution of the directions of the moments for %d rotations, Li%sF_4',N,ion{j}.name);
    title(titre)

% filename=sprintf('distrib_angles_%s.txt',ion{j}.name);
% fid=fopen(filename,'wt');
% fprintf(fid, '%s \n',colnames);
% dis_alpha=histc(alpha, angles0pi);
% dis_beta=histc(beta, angles02pi-pi);
% dis_theta=histc(theta, angles0pi-pi/2);
% dis_phi=histc(phi, angles02pi-pi);
% 
% for k=1:n
%     fprintf(fid,'%f %f %d %d %d %d \n',angles0pi(k),angles02pi(k), dis_alpha(k), dis_beta(k), dis_theta(k), dis_phi(k)); 
% end
% fclose(fid);
% end