function sout = mergeout_v1(out,opt,ind)
% Function to merge data and sort it according to optional field structure
%
% INPUT:
% out           Structure containing the data
% opt           [optional] Sort data according to this field name,
%               ie 'SerialTime', by default take first column
% ind           [optional] Indecies of "out" to merge
%               ie [1 2 3], by default take all, other


% Initiate defaults
if nargin == 1
    opt = 'SerialTime';
    ind = 1:length(out);
end
if nargin == 2 
    ind = 1:length(out);
end

% Basic check input is ok
if all(strcmpi(out(1).headers,opt))
    error('Check parameter by which to sort')
end

%% Merge data
%!!! Note that headers, instrument names, filepaths and filenames are not
% merged to retain same structure as input "out"

sout = blankstructure;

sout.headers    = out(1).headers;
sout.instrument = out(1).instrument;
sout.filepath   = out(1).filepath;
sout.filename   = out(1).filename;
    
fields = fieldnames(out(1).data);


%% Merge data from several files into one structure
for n = ind(:)'
    
    for m = 1:length(fields)
        % Concenate the data and test to see if errors are included
        % Each scan starts with timer set to 0, when concenating data, try
        % to keep timer sequential
        ff = char(fields(m));
        try
            if strcmpi(ff,'time')
                if isempty(sout.data.(ff))
                    addtime = 0;
                else
                    addtime = max(sout.data.(ff));
                end
                sout.data.(ff) = [sout.data.(ff);          out(n).data.(ff) + addtime];
            else
                sout.data.(ff) = [sout.data.(ff);          out(n).data.(ff)];
            end
        catch ME
            % If errors are not defined, replace the empty array with NaNs
            if strcmp(ff(1),'d') && ~isfield(out(n).data,ff)                
                sout.data.(ff) = [sout.data.(ff);          out(n).data.(ff(2:end)).*NaN];
            end
        end                
    end
    
    % Include list of run numbers
    sout.run        = [sout.run; out(n).run];
    
end


%% Sort the data
[sx indx] = sort(sout.data.(opt));
sout.data.(opt) = sx;

fields(strcmpi(fields,opt) == 1) = [];

for n = 1:length(fields)
    tt = sout.data.(char(fields(n)));
    sout.data.(char(fields(n))) = tt(indx,:);
end

end