% Plot the four spin configurations
temp = 0.1; % Temperature
theta = 0.0; % ab-plane angle
phi = 0.0; % c-axis angle
mion = 'Er';
hyp = 1;

filepath = ['G:\My Drive\File sharing\PhD program\Research projects\LiErF4 project\Data',...
                '\Simulations\Matlab\Susceptibilities\with Hz_I'];
filename = ['Hscan_Li',mion,'F4_',...
         sprintf('%1$3.3fK_%2$.2fDeg_%3$.1fDeg_hyp=%4$.2f',temp,theta,phi,hyp),'.mat'];
file = fullfile(filepath,filename);
load(file,'-mat','vvv','ttt','eee','fff','ion');

Jx = cell(1,4);
Jy = cell(1,4);
Jz = cell(1,4);
Jnorm = cell(1,4);

figure;
hold on
for ii = 1:4
    Jx{ii} = double.empty(4,length(fff(1,:)),0);
    Jy{ii} = double.empty(4,length(fff(1,:)),0);
    Jz{ii} = double.empty(4,length(fff(1,:)),0);
    for jj = 1:length(fff(1,:))
        Jx{ii}(:,jj,1) = ion.Js_hyp(ii,1,jj,1);
        Jy{ii}(:,jj,1) = ion.Js_hyp(ii,2,jj,1);
        Jz{ii}(:,jj,1) = ion.Js_hyp(ii,3,jj,1);
    end
    Jnorm{ii} = vecnorm([Jx{ii};Jy{ii};Jz{ii}]);
    plot(vecnorm(fff),Jnorm{ii},'o')
end
legend('Er-1','Er-2','Er-3','Er-4')
xlabel('Magnetic field (T)')
ylabel('<J>')

Jmx = double.empty(4,length(fff(1,:)),0);
Jmy = double.empty(4,length(fff(1,:)),0);
Jmz = double.empty(4,length(fff(1,:)),0);
for ii = 1:length(fff(1,:))
    Jmx(:,ii,1) = ion.Jmom_hyp(1,ii,1);
    Jmy(:,ii,1) = ion.Jmom_hyp(2,ii,1);
    Jmz(:,ii,1) = ion.Jmom_hyp(3,ii,1);
end
figure
hold on
plot(vecnorm(fff),Jmx,'-s')
plot(vecnorm(fff),Jmy,'-s')
plot(vecnorm(fff),Jmz,'-s')
legend('<Jx>','<Jx>','<Jz>')
xlabel('Magnetic field (T)')
ylabel('<J_i>')