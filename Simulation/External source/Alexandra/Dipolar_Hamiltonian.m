crystal = 'KHo3F10';
N = 3;
plot_large = 0;

if strcmp(crystal,'KHo3F10')
    a=[11.571 0      0
       0      11.571 0
       0      0      11.571];
    tau = [[0.2401, 0., 0.]; [0.7599, 0., 0.]; [0., 0.2401, 0.]; [0., 0.7599, 0.];...
           [0., 0., 0.2401]; [0., 0., 0.7599]; [0.2401, 0.5, 0.5];...
           [0.7599, 0.5, 0.5]; [0., 0.7401, 0.5]; [0., 0.2599, 0.5];...
           [0.,0.5, 0.7401]; [0., 0.5, 0.2599]; [0.7401, 0., 0.5];... 
           [0.2599, 0., 0.5]; [0.5, 0.2401, 0.5]; [0.5, 0.7599, 0.5];...
           [0.5, 0., 0.7401]; [0.5, 0., 0.2599]; [0.7401, 0.5, 0.];...
           [0.2599, 0.5, 0.]; [0.5, 0.7401, 0.]; [0.5, 0.2599, 0.];...
           [0.5, 0.5, 0.2401]; [0.5, 0.5, 0.7599]];
    mom = [[1, 0, 0]; [-1, 0, 0]; [0, -1, 0]; [0, 1, 0]; [0, 0, 1]; [0, 0, 1];...
     [1, 0, 0]; [-1, 0, 0]; [0, -1, 0]; [0, 1, 0]; [0, 0, 1]; [0, 0, 1];...
     [1, 0, 0]; [-1, 0, 0]; [0, -1, 0]; [0, 1, 0]; [0, 0, 1]; [0, 0, 1];...
     [1, 0, 0]; [-1, 0, 0]; [0, -1, 0]; [0, 1, 0]; [0, 0, 1]; [0, 0, 1]];
 sz = 100;
elseif strcmp(crystal,'LiHoF4')
a=[[5.162 0 0];[0 5.162 0];[0 0 10.70]];
    mom = [[0 0 1];[0 0 1];[0 0 1];[0 0 1]];
    tau=[[0 0 0];[0 1/2 1/4];[1/2 1/2 1/2];[1/2 0 3/4]];
sz = 600;
end

tau=tau*a;

if N>0
[x,y,z]=meshgrid(-N:N,-N:N,-N:N); % new method by HMR 2009jul - faster
hkl=[z(:) x(:) y(:)]; % use z x y to get nicer order in the list - but of no importance
r0=hkl*a;
    r1 = [];
    mom1 = [];
    for i = 1:length(r0)
        r1 = [r1;tau+repmat(r0(i,:),[size(tau,1),1])];
        mom1 = [mom1;mom];
    end
sz = sz/((2*N+1)^3);    
else
    r1 = tau;
    mom1 = mom;
end
dotmom = mom1*mom1';

momdir = zeros(size(r1,1),1);
for i = 1:size(r1,1)   
  [~,ind] = max(abs(mom1(i,:)));
   momdir(i) = ind*sign(mom1(i,ind));
end 
% close all
% figure;
% scatter3(r1(:,1),r1(:,2),r1(:,3),100,momdir,'filled')
% cb = colorbar;
%     cb.TickLabelInterpreter = 'latex';
%     cb.Label.String = '$m_i$';
%     cb.Label.Interpreter = 'latex';
% colormap(bluewhitered)
% set(gca,'ticklabelinterpreter','latex','FontSize',14)
% xlabel('$a$ (\AA)','Interpreter','latex','FontSize',14)
% ylabel('$b$ (\AA)','Interpreter','latex','FontSize',14)
% zlabel('$c$ (\AA)','Interpreter','latex','FontSize',14)

modr = zeros(size(r1,1),size(r1,1));
modmom = zeros(size(r1,1),size(r1,1));
dimom1 = zeros(size(r1,1),size(r1,1));
dimom2 = zeros(size(r1,1),size(r1,1));
dimom = zeros(size(r1,1),size(r1,1));
for i = 1:size(r1,1)
    for j = 1:size(r1,1)
        r = (r1(i,:)-r1(j,:));
        modr(i,j) = sum(r.^2);
        modmom(i,j) = (mom1(i,:)*r')*(mom1(j,:)*r');
        dimom1(i,j) = 3*(modmom(i,j)/modr(i,j)^5);
        dimom2(i,j) = -dotmom(i,j)/modr(i,j)^3;
        dimom(i,j) = dimom1(i,j)+dimom2(i,j);
    end
end

dimom(isnan(dimom) | abs(dimom)==Inf)=0;
dimom1(isnan(dimom1) | abs(dimom1)==Inf)=0;
dimom2(isnan(dimom2) | abs(dimom2)==Inf)=0;

[X,Y] = meshgrid(1:size(r1,1),1:size(r1,1));
X = reshape(X,[size(r1,1)^2 1]);
Y = reshape(Y,[size(r1,1)^2 1]);
Zmom = reshape(dotmom,[size(r1,1)^2,1]);
Zmodr = reshape(modr,[size(r1,1)^2,1]);
Zmodmom = reshape(modmom,[size(r1,1)^2,1]);
Zdimom1 = reshape(dimom1,[size(r1,1)^2,1]);
Zdimom2 = reshape(dimom2,[size(r1,1)^2,1]);
Zdimom = reshape(dimom,[size(r1,1)^2,1]);

if plot_large
    figure
    subplot(2,2,1)
    hold on
    scatter(1:size(r1,1),zeros(size(r1,1),1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(zeros(size(r1,1),1),1:size(r1,1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(0,0,sz,-3,'filled','Marker','s','MarkerEdgeColor','y')
    scatter(X,Y,sz,Zmom,'filled','Marker','s')
    if N>0
        for i = 1:size(momdir,1)/size(tau,1)
            line([-1 size(momdir,1)],ones(2,1)*(i-1)*size(tau,1),'LineStyle','-','Color','black')
            line(ones(2,1)*(i-1)*size(tau,1),[-1 size(momdir,1)],'LineStyle','-','Color','black')
        end
    end
    line([-1 size(momdir,1)],[-1 size(momdir,1)],'LineStyle','-','Color','black')
    cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.String = '$m_i \cdot m_j$';
        cb.Label.Interpreter = 'latex';
    colormap(bluewhitered)
    set(gca,'ticklabelinterpreter','latex','FontSize',14)
    xlabel('Ion Number','Interpreter','latex','FontSize',14)
    ylabel('Ion Number','Interpreter','latex','FontSize',14)
    xlim([-1 size(r1,1)+1])
    ylim([-1 size(r1,1)+1])

    % figure;
    subplot(2,2,2)
    hold on
    scatter(1:size(r1,1),zeros(size(r1,1),1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(zeros(size(r1,1),1),1:size(r1,1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(0,0,sz,-3,'filled','Marker','s','MarkerEdgeColor','y')
    scatter(X,Y,sz,Zdimom2/max(abs(Zdimom2)),'filled','Marker','s')
    if N>0
        for i = 1:size(momdir,1)/size(tau,1)
            line([-1 size(momdir,1)],ones(2,1)*(i-1)*size(tau,1),'LineStyle','-','Color','black')
            line(ones(2,1)*(i-1)*size(tau,1),[-1 size(momdir,1)],'LineStyle','-','Color','black')
        end
    end
    line([-1 size(momdir,1)],[-1 size(momdir,1)],'LineStyle','-','Color','black')
    cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.String = '$\frac{1}{|r_{ij}|^3}(m_i \cdot m_j)$';
        cb.Label.Interpreter = 'latex';
    colormap(bluewhitered)
    set(gca,'ticklabelinterpreter','latex','FontSize',14)
    xlabel('Ion Number','Interpreter','latex','FontSize',14)
    ylabel('Ion Number','Interpreter','latex','FontSize',14)
    xlim([-1 size(r1,1)+1])
    ylim([-1 size(r1,1)+1])

    % figure;
    subplot(2,2,3)
    hold on
    scatter(1:size(r1,1),zeros(size(r1,1),1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(zeros(size(r1,1),1),1:size(r1,1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(0,0,sz,-3,'filled','Marker','s','MarkerEdgeColor','y')
    scatter(X,Y,sz,Zdimom1/max(abs(Zdimom1)),'filled','Marker','s')
    if N>0
        for i = 1:size(momdir,1)/size(tau,1)
            line([-1 size(momdir,1)],ones(2,1)*(i-1)*size(tau,1),'LineStyle','-','Color','black')
            line(ones(2,1)*(i-1)*size(tau,1),[-1 size(momdir,1)],'LineStyle','-','Color','black')
        end
    end
    line([-1 size(momdir,1)],[-1 size(momdir,1)],'LineStyle','-','Color','black')
    cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.String = '$\frac{3}{|r_{ij}|^5}(m_i \cdot r_{ij}) * (m_j \cdot r_{ij})$';
        cb.Label.Interpreter = 'latex';
    colormap(bluewhitered)
    set(gca,'ticklabelinterpreter','latex','FontSize',14)
    xlabel('Ion Number','Interpreter','latex','FontSize',14)
    ylabel('Ion Number','Interpreter','latex','FontSize',14)
    xlim([-1 size(r1,1)+1])
    ylim([-1 size(r1,1)+1])

% figure;
    subplot(2,2,4)
    hold on
    scatter(1:size(r1,1),zeros(size(r1,1),1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(zeros(size(r1,1),1),1:size(r1,1),sz,momdir,'filled','Marker','s','MarkerEdgeColor','black')
    scatter(0,0,sz,-3,'filled','Marker','s','MarkerEdgeColor','y')
    scatter(X,Y,sz,Zdimom/max(abs(Zdimom)),'filled','Marker','s')
    if N>0
        for i = 1:size(momdir,1)/size(tau,1)
            line([-1 size(momdir,1)],ones(2,1)*(i-1)*size(tau,1),'LineStyle','-','Color','black')
            line(ones(2,1)*(i-1)*size(tau,1),[-1 size(momdir,1)],'LineStyle','-','Color','black')
        end
    end
    line([-1 size(momdir,1)],[-1 size(momdir,1)],'LineStyle','-','Color','black')
    cb = colorbar;
        cb.TickLabelInterpreter = 'latex';
        cb.Label.String = '$\frac{3}{|r_{ij}|^5}(m_i \cdot r_{ij}) * (m_j \cdot r_{ij})-\frac{1}{|r_{ij}|^3}(m_i \cdot m_j)$';
        cb.Label.Interpreter = 'latex';
    colormap(bluewhitered)
    set(gca,'ticklabelinterpreter','latex','FontSize',14)
    xlabel('Ion Number','Interpreter','latex','FontSize',14)
    ylabel('Ion Number','Interpreter','latex','FontSize',14)
    xlim([-1 size(r1,1)+1])
    ylim([-1 size(r1,1)+1])
end

rr = sum(r1.^2,2);
dimom_sum = sum(dimom,1);
ind_sph = find(rr<=max(rr)*0.1);

[x, y, z] = sphere;
x = x * max(rr)*0.1;
y = y * max(rr)*0.1;
z = z * max(rr)*0.1;



figure;
hold on
        plot3([0 1]*max(a(1,:)),[0 0]*max(a(2,:)),[0 0]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 0]*max(a(1,:)),[0 1]*max(a(2,:)),[0 0]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 0]*max(a(1,:)),[0 0]*max(a(2,:)),[0 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([1 1]*max(a(1,:)),[0 0]*max(a(2,:)),[0 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 0]*max(a(1,:)),[1 1]*max(a(2,:)),[0 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([1 1]*max(a(1,:)),[1 1]*max(a(2,:)),[0 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([1 1]*max(a(1,:)),[0 1]*max(a(2,:)),[0 0]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 0]*max(a(1,:)),[0 1]*max(a(2,:)),[1 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([1 1]*max(a(1,:)),[0 1]*max(a(2,:)),[1 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 1]*max(a(1,:)),[1 1]*max(a(2,:)),[0 0]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 1]*max(a(1,:)),[0 0]*max(a(2,:)),[1 1]*max(a(3,:)),'-black','LineWidth',2)
        plot3([0 1]*max(a(1,:)),[1 1]*max(a(2,:)),[1 1]*max(a(3,:)),'-black','LineWidth',2)
scatter3(r1(ind_sph,1),r1(ind_sph,2),r1(ind_sph,3),100,dimom_sum(ind_sph),'filled')
cb = colorbar;
    cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$H_{dipole}$';
    cb.Label.Interpreter = 'latex';
colormap(bluewhitered)
pbaspect([max(a(1,:)),max(a(2,:)),max(a(3,:))])
set(gca,'ticklabelinterpreter','latex','FontSize',14)
xlabel('$a$ (\AA)','Interpreter','latex','FontSize',14)
ylabel('$b$ (\AA)','Interpreter','latex','FontSize',14)
zlabel('$c$ (\AA)','Interpreter','latex','FontSize',14)

