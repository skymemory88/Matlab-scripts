function tout = sub_v1(out1,field1,out2,field2)
% Function to subtract two datasets from one another as out1 - out2
% Note that if the objects are different length, will not subtract
% Returned object out will contain the subtraction in "field1" and "field2"

if length(out1) ~= length(out2) || length(out1) > 1 || length(out2) > 1
    error('Input data must be 1D')    
end

% Extract data from structures
[y1 dy1] = extractvalues(out1,field1);
[y2 dy2] = extractvalues(out2,field2);

% Check that arrays are the same length before trying subtraction
if length(y1) ~= length(y2)
    error('Objects different lengths')    
end

% Subtract one dataset from another and calculate the errors
y12 = y1 - y2;
dy12 = sqrt(dy1.^2 + dy2.^2);

% Load a blank file structure
tout = blankstructure;
fnames = fieldnames(tout.data);
Nf = length(fnames);

for n = 1:Nf/2
    tout.data.(fnames{n}) = mean([out1.data.(fnames{n}) out2.data.(fnames{n})],2);
    if ~isempty(out1.data.([ 'd' fnames{n}])) && ~isempty(out1.data.([ 'd' fnames{n}]))
        tout.data.([ 'd' fnames{n}]) = sqrt(out1.data.([ 'd' fnames{n}]).^2 + out2.data.([ 'd' fnames{n}]).^2);    
    end    
end

% Load the subtraction into the data file
tout.data.(field1) = y12;
tout.data.(field2) = y12;

tout.data.(['d' field1]) = dy12;
tout.data.(['d' field2]) = dy12;

% figure, plot(y12)
end