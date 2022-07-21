% constructing phase diagram
% temp = [0.05 [0.1:0.05:0.4] 0.45:0.05:1.65 1.675 1.7]; % temperature points
temp = [0.05 0.1:0.1:1.6 1.65:0.05:1.8];

Options.ion = 'Ho'; % LiReF4: Re ion selection
    elem_idx = find(strcmp(Options.ion,[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}]));
Options.nZee = false;
Options.hyp = 1;

theta = 0.0;
phi = 4.0;

lin = ['-','--','-.',':']; % line style pool

% figure
box on
hold on

if Options.nZee == true
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',Options.ion,...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I'];
else
    Options.location = ['G:\My Drive\File sharing\PhD program\Research projects\Li',Options.ion,...
        'F4 project\Data\Simulations\Matlab\Susceptibilities\without Hz_I\R=0.785'];
end

bound = double.empty(length(temp), 2, 0); % container for the phase boundary
for ii = 1:length(temp)
    lname=['Hscan_Li',Options.ion,'F4_', sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', temp(ii), theta, phi, Options.hyp),'.mat'];
    file = fullfile(Options.location, lname);
    load(file,'-mat', 'fff', 'ion');
%     load(file,'-mat','vvv','ttt','eee','fff','ion');
%     eee = eee- min(eee,[],2); % Normalize the energy amplitude to the lowest eigen-energy
    field = vecnorm(fff,1); % take the vectorial norm of the external magnetic field
    Js = squeeze(ion.Jmom_hyp(:,:,elem_idx));
    [~,idx] = find(Js(3,:) <= 1e-3, 1, 'first');
    bound(ii,1,1) = temp(ii);
    bound(ii,2,1) = field(idx);
end
% plot(bound(:,1), bound(:,2),'o');
phase = plot(bound(:,1), bound(:,2), 'LineStyleMode', 'manual', 'Color', 'r', 'LineStyle', lin(1));
xlabel('Temperature (K)')
ylabel('Magnetic field (T)')