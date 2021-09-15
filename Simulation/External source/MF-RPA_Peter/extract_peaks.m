function s = extract_peaks(s, xvec, yvec, zvec, dim, show, ztol, ytol, opt_sort)
% Note: xvec, yvec, zvec are figure axes

extpeaks.pks = zeros(length(s.(xvec)),20).*NaN;
extpeaks.loc = zeros(length(s.(xvec)),20).*NaN;
extpeaks.wid = zeros(length(s.(xvec)),20).*NaN;
extpeaks.amp = zeros(length(s.(xvec)),20).*NaN;

ycut = s.(yvec);

try
    ind_ycut = ycut>ytol(1) & ycut<ytol(2);
    ycut = ycut(ind_ycut);
end

for n=1:length(s.(xvec))
    tt = squeeze(s.(zvec));
    if dim == 1
        zcut = tt(n,:);
    else
        zcut = tt(:,n);
    end
    
    try
        zcut = zcut(ind_ycut);
    end
    
    [pks,loc,wid,amp] = findpeaks(zcut,ycut,'MinPeakProminence',ztol);
    N = length(pks);
    extpeaks.pks(n,1:N) = pks;
    extpeaks.loc(n,1:N) = loc;
    extpeaks.wid(n,1:N) = wid;
    extpeaks.amp(n,1:N) = amp;
    
    try
        ind_weak = pks < ztol;
        pks(ind_weak) = [];
        loc(ind_weak) = [];
        wid(ind_weak) = [];
        amp(ind_weak) = [];
    end
    
    ind_inf = isinf(pks);
    pks(ind_inf) = [];
    loc(ind_inf) = [];
    wid(ind_inf) = [];
    amp(ind_inf) = [];
    N = length(pks);
    
    try
        if ~opt_sort
            ind = 1:N;
        else
            [dum,ind] = sort(pks,'descend');
        end
    end
        
    extpeaks.spks(n,1:N) = pks(ind);
    extpeaks.sloc(n,1:N) = loc(ind);
    extpeaks.swid(n,1:N) = wid(ind);
    extpeaks.samp(n,1:N) = amp(ind);
    extpeaks.iint(n) = trapz(ycut,zcut);
    extpeaks.xpar(n) = s.(xvec)(n);
    
%     try
    if show == true
        figure(999)
        clf
        plot(ycut,zcut)
        hold on
%         plot(extpeaks.sloc(n,:),extpeaks.spks(n,:),'or')
%         text(extpeaks.sloc(n,:),extpeaks.spks(n,:),num2str([1:length(extpeaks.sloc(n,:))]'))
        plot(loc(ind),pks(ind),'or')
        text(loc(ind),pks(ind),num2str([1:length(ind)]'))
        title(num2str([n,length(s.(xvec))],'%d / %d'))
        
        pause(0.5)
        %keyboard
    end
%     end
    
end

extpeaks.xvec = xvec;
extpeaks.yvec = yvec;
extpeaks.zvec = zvec;
extpeaks.dim = dim;

s.extpeaks = extpeaks;

end