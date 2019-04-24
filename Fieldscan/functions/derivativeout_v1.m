function out = derivativeout_v1(out,varx,vary)
% Function to take the numerical derivative of vary with respect to varx
%
% INPUT:
% out:                              structure containing data
% varx:                             differntiation wrt to this parameter
% vary:                             function to be differentiated
%
% OUTPUT:
% out.data.(['diff' varx]):        derivatives only structure out containing x axis
% out.data.(['diff' vary]):        derivatives only structure out containing dy/dx



for n = 1:length(out)
    x = out(n).data.(varx);
    y = out(n).data.(vary);
    
    
    
    out(n).data.(['diff_' varx]) = x;
    out(n).data.(['diff_' vary]) = diff(y)./diff(x);
    
end

end