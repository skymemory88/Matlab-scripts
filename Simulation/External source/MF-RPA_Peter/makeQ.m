function q = makeQ(h,k,l)

% Input:
% h = 0;
% k = -0.2:0.05:0.2;
% l = -1;

Nh = length(h);
Nk = length(k);
Nl = length(l);

maxN = max([Nh Nk Nl]);

if Nh == 1
    qh = h*ones(maxN,1);
else
    qh = h(:);
end

if Nk == 1
    qk = k*ones(maxN,1);
else
    qk = k(:);
end

if Nl == 1
    ql = l*ones(maxN,1);
else
    ql = l(:);
end


q = [qh qk ql];


end