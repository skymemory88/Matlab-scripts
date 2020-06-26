function [C_fdt_fwd, C_fdt_bwd, elems] = specHt(E, params)
    temps = size(E,1);
    kB = 8.617E-2; % Boltzmann constant [meV/K]
    totN = floor(size(E,2)/params.N_Er);
    C_fdt_fwd = double.empty(temps,totN-1,0);
    C_fdt_bwd = double.empty(temps,totN-1,0);
    elems = double.empty(totN-1,0);
    for ii = 1:temps
        for kk = 1:totN-1
        C_fdt_fwd(ii,kk,1) = (mean(E(ii,1:(kk+1)*params.N_Er).^2) - mean(E(ii,1:(kk+1)*params.N_Er))^2)/params.temp(ii)/kB;
        C_fdt_bwd(ii,totN-kk,1) = (mean(E(ii,kk*params.N_Er:end).^2) - mean(E(ii,kk*params.N_Er:end))^2)/params.temp(ii)/kB;
        elems(kk,1) = (kk+1)*params.N_Er;
        end
    end
end