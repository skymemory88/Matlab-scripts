function pos = lattice(uVec, bVec, params)
    % Construct a lattice with a finite size and custom cutoff shape
    % Input:
    %   bVec - n x 3 matrix containing the basis vectors as rows
    %   uVec - n x 3 matrix containing the unit vectors as rows
    %   dims - Vector [nx, ny, nz] specifying the number of unit cells along each axis
    %   cutR - Cutoff radius for the spherical cutoff
    % Output:
    %   pos - n x 3 matrix containing the x, y, z coordinates of lattice points
    
    dims = params.dims; % lattice dimension
    cutR = params.domRng/2 * min(diag(uVec),[],'all'); % lattice cutoff radius
    maskType = params.coord; % shaping mask type ('spherical', 'cubic', 'ellipsoidal')

    % Generate unit cell indices centered around zero
    [X, Y, Z] = ndgrid(-floor(dims(1)/2):ceil(dims(1)/2)-1, ...
                       -floor(dims(2)/2):ceil(dims(2)/2)-1, ...
                       -floor(dims(3)/2):ceil(dims(3)/2)-1);
    
    % Reshape to vectors
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    % Generate base positions for all unit cells
    basePos = [X, Y, Z] * uVec;
    
    % correct the unit vector as they are defined by ratio instead of absolute length
    bVec = bVec * uVec;
    
    % Generate all lattice points based on basis vectors
    allPos = bsxfun(@plus, permute(basePos, [1, 3, 2]), permute(bVec, [3, 1, 2]));
    allPos = reshape(allPos, [], 3);
    
    % Apply custom cutoff shape
    cutoffMask = arrayfun(@(x, y, z) cutFunc([x, y, z], maskType, cutR), allPos(:, 1), allPos(:, 2), allPos(:, 3));
    
    % Final positions
    pos = allPos(cutoffMask, :);
end

function isInside = cutFunc(position, cutType, varargin)
    % cutFunc: A function to apply different types of cutoffs
    % position: The position vector to be checked
    % cutType: The type of cutoff ('spherical', 'cubic', 'ellipsoidal')
    % varargin: Additional parameters required for the cutoff
    
    switch cutType
        case 'spherical'
            if numel(varargin) < 1
                error('Radius is required for spherical cutoff.');
            end
            radius = varargin{1};
            isInside = norm(position) <= radius;
            
        case 'cubic'
            if numel(varargin) < 1
                error('Side length is required for cubic cutoff.');
            end
            sideLength = varargin{1};
            isInside = all(abs(position) <= sideLength / 2);
            
        case 'ellipsoidal'
            if numel(varargin) < 3
                error('Semi-axes a, b, and c are required for ellipsoidal cutoff.');
            end
            a = varargin{1};
            b = varargin{2};
            c = varargin{3};
            x = position(1);
            y = position(2);
            z = position(3);
            isInside = (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1;
            
        otherwise
            error('Invalid cutType. Choose from "spherical", "cubic", or "ellipsoidal".');
    end
end
