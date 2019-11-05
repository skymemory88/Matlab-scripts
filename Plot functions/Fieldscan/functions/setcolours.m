function col = setcolours(perc,type)
% Function to generate a colour scheme for plots
%
% INPUT:
% perc                  fraction of the colour vector to use
% type                  different types of plot colour schemes
%
% OUTPUT:
% col                   1x3 vector containing rgb colour information

try
    cmap = colormap(type);
catch ME
    if ~isstrprop(type,'digit')
        type = 'jet';
        cmap = colormap(type);
    end
end

switch lower(type)
        
    case 'rand'
        col = rand(1,3);    
    case '1'
        col = [0.7 0.2 0.2];
    case '2'
        col = [0.2 0.7 0.2];
    case '3'
        col = [0.2 0.2 0.7];
    case '4'
        col = [0.7 0.2 0.7];
    otherwise
        % Load the specified colour map
        N = length(cmap);
        col = interp1(cmap,perc*(N-1)+1 );

end

end
