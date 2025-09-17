classdef GPUSpinState < handle
    properties
        % Persistent GPU arrays (shared across all temperatures)
        d_pos           % positions on GPU
        d_nList         % neighbor lists on GPU (cell array)
        d_nVecs         % neighbor vectors on GPU
        d_nDists        % neighbor distances on GPU
        d_gfac          % dipole prefactor
        
        % Per-temperature GPU arrays
        d_eSpin         % spin configurations [temp, field]
        
        % Pair data for optimized dipSum calculation
        pairData        % struct containing pre-built pair information
        
        % System parameters
        N               % number of spins
        nTemps          % number of temperatures
        nFields         % number of fields
        initialized     % initialization flag
        use_double      % precision flag
    end
    
    methods
        function obj = GPUSpinState(params, nTemps, nFields, use_double)
            obj.N = size(params.pos, 1);
            obj.nTemps = nTemps;
            obj.nFields = nFields;
            obj.initialized = false;
            obj.pairData = [];  % Initialize as empty
            if nargin < 4
                obj.use_double = true;  % Default to double precision for accuracy
            else
                obj.use_double = use_double;
            end
        end
        
        function initialize(obj, params, gfac)
            if ~obj.initialized
                if obj.use_double
                    precision = 'double';
                else
                    precision = 'single';
                end
                
                % Transfer persistent data with consistent precision
                obj.d_pos = gpuArray(cast(params.pos * 1e-10, precision));
                obj.d_gfac = gpuArray(cast(gfac, precision));
                
                % Convert neighbor lists to GPU format
                obj.d_nList = cell(obj.N, 1);
                obj.d_nVecs = cell(obj.N, 1);
                obj.d_nDists = cell(obj.N, 1);
                
                for i = 1:obj.N
                    if ~isempty(params.nList{i})
                        obj.d_nList{i} = gpuArray(int32(params.nList{i}));
                        obj.d_nVecs{i} = gpuArray(cast(params.nVecs{i}, precision));
                        obj.d_nDists{i} = gpuArray(cast(params.nDists{i}, precision));
                    end
                end
                
                % Allocate per-temperature arrays
                obj.d_eSpin = cell(obj.nTemps, obj.nFields);
                for t = 1:obj.nTemps
                    for f = 1:obj.nFields
                        obj.d_eSpin{t,f} = gpuArray.zeros(obj.N, 3, precision);
                    end
                end
                
                obj.initialized = true;
                fprintf('GPU state initialized: %d spins, %d temps, %d fields (%s precision)\n', ...
                    obj.N, obj.nTemps, obj.nFields, precision);
            end
        end
        
        function updateSpins(obj, eSpin, tempIdx, fieldIdx)
            if obj.use_double
                obj.d_eSpin{tempIdx, fieldIdx} = gpuArray(double(eSpin));
            else
                obj.d_eSpin{tempIdx, fieldIdx} = gpuArray(single(eSpin));
            end
        end
        
        function eSpin = getSpins(obj, tempIdx, fieldIdx)
            eSpin = gather(obj.d_eSpin{tempIdx, fieldIdx});
        end
    end
end