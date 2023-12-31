% constructing phase diagram
% temp = [0.05 0.1:0.1:1.4]; % temperature points
temp = [0 0.05 0.1:0.1:1.6 1.65:0.05:2.0];

Options.ion = 'Ho'; % LiReF4: Re ion selection
    elem_idx = find(strcmp(Options.ion,[{'Er'};{'Ho'};{'Yb'};{'Tm'};{'Gd'};{'Y'}]));
Options.nZee = true;
Options.hyp = 1.00;
Options.sim = true;
    theta = 0.0;
    phi = 12.0;
Options.exp = false;
    Options.sample = 239;
%     Options.sample = 200;

lin = ["-" "--", "-.", ":", "-*"]; % line style pool
mkr = ["o", "s", "x", "*", "^"]; % marker style pool
col = ["k", "r", "b", "m", "g"]; % built-in color pool
mc = {[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]; ...
    [0.4940, 0.1840, 0.5560]; [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330]}; % marker color pool

% figure
box on
hold on
if Options.sim == true
    if Options.nZee == true
        sim_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\Research projects\Li',...
            Options.ion, 'F4 project\Data\Simulations\Matlab\Susceptibilities\with Hz_I\phase'];
    else
        sim_loc = ['G:\.shortcut-targets-by-id\1CapZB_um4grXCRbK6t_9FxyzYQn8ecQE\File sharing\PhD program\Research projects\Li',...
            Options.ion, 'F4 project\Data\Simulations\Matlab\Susceptibilities\without Hz_I\phase'];
    end

    bound = double.empty(length(temp), 2, 0); % container for the phase boundary
    for ii = 1:length(temp)
        lname=['Hscan_Li',Options.ion,'F4_', sprintf('%1$3.3fK_%2$.2fDg_%3$.1fDg_hp=%4$.2f', temp(ii), theta, phi, Options.hyp),'.mat'];
        file = fullfile(sim_loc, lname);
        load(file,'-mat', 'fff', 'ion');
%         load(file,'-mat','vvv','ttt','eee','fff','ion');
%         eee = eee- min(eee,[],2); % Normalize the energy amplitude to the lowest eigen-energy
        field = vecnorm(fff,2,1); % take the vectorial norm of the external magnetic field
        Js = squeeze(ion.Jmom_hyp(:,:,elem_idx));
        if theta == 0
            [~,idx] = find(Js(3,:) <= 1e-4, 1, 'first');
        else
            disp('Use non-vaishing order parameter threshold\n')
            [~,idx] = find(Js(3,:) <= 1.4, 1, 'first'); % only for when there is a longitudinal field
        end
        bound(ii,1,1) = temp(ii);
        bound(ii,2,1) = field(idx);
    end
    sim_phase = plot(bound(:,1), bound(:,2), 'LineStyle', lin(1), 'Color', col(2));
    xlabel('Temperature (K)')
    ylabel('Magnetic Field (T)')
end

if Options.exp == true
    exp_loc = ['G:\My Drive\File sharing\PhD program\Research projects\LiHoF4 project\Data\Experiment\LiHoF4\SC',...
        num2str(Options.sample), '\'];
    dat_lst = readtable(fullfile(exp_loc, sprintf("SC%u_off_list.xlsx", Options.sample)));
    folder = table2array(dat_lst(:,1));
    efiles = table2array(dat_lst(:,2));
    etemps = table2array(dat_lst(:,3));
    B0 = double.empty(length(etemps),0);
    for ii = 1:length(etemps)
        exp_file = fullfile([exp_loc, folder{ii},'\'], efiles{ii}); % locate the experimental data
        load([exp_file, '_interp'], '-mat', 'analysis', 'continu_var');
        [~,idx] = min(analysis.w0);
        B0(ii,1) = analysis.Hs(idx);
    end
    exp_phase = scatter(etemps, B0, 20, 'Marker', mkr(1), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
end