function tout = times_v1(out1,field1,multfac)
% Function to multiply field1 of structure 1 by factor multfac
% Note that if the objects are different length, will not multiply
% Returned object out will contain the subtraction in "field1"
% Modified: PB 28.03.2014

% Extract data from structure
[y1 dy1] = extractvalues(out1,field1);

% Subtract one dataset from another and calculate the errors
y12 = y1.*multfac;
dy12 = dy1.*multfac;

tout = out1;
for n=1:length(out1)
    % Load the subtraction into the data file
    tout(n).data.(field1) = y12;
    tout(n).data.(['d' field1]) = dy12;
end
    

end