function para_diag_nohyp
clear all

Hcf=cf_Ho;
I=3.5;
J=8;
A=0.003361;
gLande=1.25;

%Initiate J operators
Jz=diag(J:-1:-J);
Jp=diag(sqrt((J-[(J-1):-1:-J]).*(J+1+[(J-1):-1:-J])),1);
Jm=Jp';
Jx=(Jp+Jm)/2;
Jy=(Jp-Jm)/2i;

b=0:0.01:5;
B=[b;zeros(size(b));zeros(size(b))];

for h=1:size(B,2)
    Hzeeman=(-gLande*0.05788)*(B(1,h)*Jx+B(2,h)*Jy+B(3,h)*Jz);
    Ham=Hcf+Hzeeman;
    [v,d]=eig(Ham);
%     e(h,:)=diag(d)-min(diag(d));
    e(h,:)=diag(d);
    one(h,:)=e(h,:)-e(h,1);
    two(h,:)=e(h,:)-e(h,2);
    three(h,:)=e(h,:)-e(h,3);
    jz(h,:)=diag(v'*Jz*v);
    jx(h,:)=diag(v'*Jx*v);
    jy(h,:)=diag(v'*Jy*v);
end
hplanck = 4.13566733e-15; % eV*s  is Planck's constant
% e = e / 1e3 / hplanck / 1e9 + 5100; % GHz
% e = e / 1e3 / hplanck / 1e9; % GHz

one = abs( one / 1e3 / hplanck / 1e9 ); % GHz
two = abs( two / 1e3 / hplanck / 1e9 ); % GHz
three = abs( three / 1e3 / hplanck / 1e9 ); % GHz
    
% figure
% plot(B(1,:),jx(:,1:16))
% figure
% plot(B(1,:),jy(:,1:16))
% figure
% plot(B(1,:),jz(:,1:16))

hfig1 = figure(4349);
clf
ttitle='para diag nohyp - Hz - Energies -E(first)';
% plot(B(3,:),e(:,1:4));
hold on
% plot(B(3,:),one(:,1:4));
% plot(B(3,:),two(:,1:4));
% plot(B(3,:),three(:,1:4));
plot(B(1,:),e(:,:));
set(hfig1,'position',[50 50 950 950])
xlabel('Hz [T]','FontSize',20)
ylabel('Frequency [GHz]','FontSize',20)
title(ttitle,'FontSize',20)
saveplots(hfig1,ttitle)

end


function saveplots(hfig,figname)
% if strcmpi(input(['Save plot ' num2str(hfig) ' {' figname '}? Y/[N]: '],'s'),'y')
curdir = cd;
cd('C:\Users\ikovacev\Documents\Projects\LiHoF4_EPR_MF_calculations_and_literature\Ivan')

saveas(figure(hfig),[figname '.fig'],'fig');
print(figure(hfig),[figname  '.jpg'],'-djpeg','-r600');
print2eps(figname,hfig)
[result,msg] = eps2xxx([figname '.eps'],{'jpeg','pdf'});

disp(['Figure ' figname ' saved to '])
disp(cd)
cd(curdir)
end
