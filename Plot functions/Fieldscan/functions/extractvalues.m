function [x dx] = extractvalues(out,varx)
% Function to extract values of parameter varx and output them as arrays x
% and dx containing data values and uncertainties, respectively.

x = [];
dx = [];

for n = 1:length(out)
    x = [x; out(n).data.(varx)];
    
    if isfield(out(n).data,['d' varx])
        dx = [dx; out(n).data.(['d' varx])];
    else
        dx = 0.*x;
    end
end


end