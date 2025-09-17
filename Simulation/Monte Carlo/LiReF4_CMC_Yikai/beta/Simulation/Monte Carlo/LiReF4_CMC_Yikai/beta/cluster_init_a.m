function params = cluster_init_a(params, const, opts)
%LB_PRECOMPUTE_TABLES_A  Build Luijten–Blöte tables and store in params.LB
% Uses params.pos (N×3) in Angstrom. If you have PBC, also set params.box = [Lx Ly Lz] in Angstrom.
% Internally converts to meters for consistent units with gfac.
% const.gfac (D) is only used to sanity-check; tables are unit-scaled (no beta/D).

    if ~isfield(params,'pos'); error('params.pos (N×3) is required'); end
    r = params.pos * 1e-10;  % Convert from Angstrom to meters for consistent units
    N = size(r,1);

    % Options (all optional)
    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'dropTies'),    opts.dropTies = true;   end
    if ~isfield(opts,'keyDigits'),   opts.keyDigits = 12;    end  % hash rounding
    if ~isfield(opts,'usePBC'),      opts.usePBC = isfield(params,'box'); end

    if opts.usePBC
        if ~isfield(params,'box'), error('Set params.box = [Lx Ly Lz] for PBC'); end
        box = params.box(:).' * 1e-10;  % Convert box dimensions from Angstrom to meters
    else
        box = [0 0 0]; % no wrapping
    end

    % --- Build relative offsets from reference site 1 to all others (wrapped if PBC) ---
    Off = zeros(N-1, 3);      % wrapped physical offsets Δr [meters]
    invr3 = zeros(N-1, 1);    % 1 / |Δr|^3 [1/meters^3]
    keep  = true(N-1, 1);     % drop ambiguous L/2 ties if requested
    k = 0;
    r1 = r(1,:);
    for j = 2:N
        k = k + 1;
        dr = r(j,:) - r1;
        dr = mi_wrap(dr, box);   % minimal image if PBC, passthrough otherwise
        if opts.dropTies && any(abs(abs(dr) - box/2) < 1e-12 & box>0)
            keep(k) = false;
        end
        Off(k,:)   = dr;
        rn = norm(dr);
        if rn == 0, invr3(k) = 0; keep(k) = false;
        else,       invr3(k) = 1 / (rn^3);
        end
    end

    % Keep only valid offsets
    Off   = Off(keep,:);
    invr3 = invr3(keep);

    % Fixed order for LB lookup; sorting by descending 1/r^3 is a good default
    [invr3_ord, ord] = sort(invr3, 'descend');
    Off_ord = Off(ord,:);

    % Unit-scaled cumulative sum: S_unit(n) = sum_{l<=n} (1/r_l^3)
    S_unit = cumsum(invr3_ord);

    % Build a robust position→index map (wrapped positions → site id)
    posmap = build_pos_map(r, box, opts.keyDigits);

    % Pack in params.LB
    LB = struct();
    LB.N          = N;
    LB.box        = box;
    LB.usePBC     = opts.usePBC;
    LB.dropTies   = opts.dropTies;
    LB.keyDigits  = opts.keyDigits;
    LB.Offsets    = Off_ord;       % (N-1)×3 [meters]
    LB.inv_r3     = invr3_ord;     % (N-1)×1 [1/meters^3]  
    LB.S_unit     = S_unit;        % (N-1)×1 [1/meters^3] (cumsum of inv_r3)
    LB.posmap     = posmap;        % containers.Map (hash of wrapped pos → index)

    params.LB = LB;

    % (Optional) small report
    if isfield(const,'gfac')
        fprintf('LB tables built: N=%d, offsets=%d, gfac=%.6g, PBC=%d\n', ...
            N, numel(invr3_ord), const.gfac, LB.usePBC);
    else
        fprintf('LB tables built: N=%d, offsets=%d, PBC=%d\n', ...
            N, numel(invr3_ord), LB.usePBC);
    end
end

% -------- helpers ----------
function rij = mi_wrap(dr, box)
% minimal-image wrap if box>0, otherwise passthrough
    rij = dr;
    for a = 1:3
        L = box(a);
        if L > 0
            rij(a) = rij(a) - L*round(rij(a)/L);
        end
    end
end

function posmap = build_pos_map(r, box, keyDigits)
    N = size(r,1);
    keys = cell(N,1);
    for i = 1:N
        riw = wrap_pos(r(i,:), box);
        keys{i} = hash_pos(riw, keyDigits);
    end
    posmap = containers.Map(keys, num2cell(1:N));
end

function rp = wrap_pos(r, box)
    rp = r;
    for a = 1:3
        L = box(a);
        if L > 0
            rp(a) = r(a) - L*floor((r(a)/L) + 0.5); % center-based wrap
        end
    end
end

function key = hash_pos(ri, keyDigits)
    fmt = sprintf('%%.%df,%%.%df,%%.%df', keyDigits, keyDigits, keyDigits);
    key = sprintf(fmt, ri(1), ri(2), ri(3));
end
