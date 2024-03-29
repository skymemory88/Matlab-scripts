function save_fit(mfields,freq,S11,analysis,H1,H2,dType)

[~,id1] = min(abs(mfields(1,:) - H1));
[~,id2] = min(abs(mfields(1,:) - H2));
switch dType
    case 'exp'
        savename = sprintf('%1$s_%2$.3fK_%3$.2f-%4$.2fT.mat',analysis.name, analysis.temp, analysis.field_l, analysis.field_h);
    case 'sim'
        savename = sprintf('S11_LiHoF4_%1$.3fK_%2$.1f-%3$.1fT_analysis.mat',analysis.temp, analysis.field_l, analysis.field_h);
end

savefile = fullfile(analysis.location,savename);

attn = analysis.attn(id1:id2);
temp = analysis.temp;
H0 = mfields(1,id1:id2)';

w0 = analysis.w0(id1:id2);
w0_ci = analysis.w0_ci(id1:id2,:);
% wc = analysis.wc(id1:id2);
% wc_ci = analysis.wc_ci(id1:id2,:);
wr = analysis.wr(id1:id2);

gc = analysis.Gc(id1:id2);
gc_ci = analysis.Gc_ci(id1:id2,:);
gc_err = abs(analysis.gma_ci(id1:id2,1)-analysis.gma_ci(id1:id2,2));

gma = analysis.gma(id1:id2);
gma_ci = analysis.gma_ci(id1:id2,:);
gma_err = abs(analysis.gma_ci(id1:id2,1)-analysis.gma_ci(id1:id2,2));

if exist('gma0','var')
    gma0 = analysis.gamma0;
    save(savefile,'mfields','freq','S11','H0','w0','w0_ci','wr','gc','gc_ci','gc_err','gma','gma_err',...
        'gma_ci','gma_err','attn','gma0','temp');
else
    save(savefile,'mfields','freq','S11','H0','w0','w0_ci','wr','gc','gc_ci','gc_err','gma','gma_err',...
        'gma_ci','gma_err','attn','temp');
end