function out = filterout_v1(out,vary,...      % necessary parameters
                thresh,varx,typ)           % optional parameters
%            
% Function to filter the data from anomalous points using derivatives of
% varx and vary based on some threshold defined by thresh
%
% INPUT:
% out:          Input structure containing fields defined in blankstructure 
%               function
% vary:         Define by which which field to filter the data, ie 'DCField1'
% thresh:       [optional] Derivative threshold, by default take 95% of
%               data and disgard 5% with largest deviation
% varx:         [optional] By default differentiate wrt "Time"
% typ:          [optional] Choose between 'abs' or 'dydx' methods for
%               filtering
%
% NB The order of the optional input does not matter for this function
%
% OUTPUT:
% fout          Return structure defined in blankstructure function



%% Initialise all the parameters
if nargin < 3
    thresh = [];
end
if nargin < 4
    varx = [];
end
if nargin < 5
    typ = [];
end
sthresh = [];
svarx = [];
styp = [];

if isnumeric(thresh) && ~isempty(thresh)
    sthresh = thresh;
else
    if strcmpi(thresh,'abs') || strcmpi(thresh,'dydx') || strcmpi(thresh,'polyfit') || strcmpi(thresh,'conv')
        styp = thresh;
    else
        svarx = thresh;
    end
end
if isnumeric(varx) && ~isempty(varx)
    sthresh = varx;
else
    if strcmpi(varx,'abs') || strcmpi(varx,'dydx') || strcmpi(varx,'polyfit') || strcmpi(varx,'conv')
        styp = varx;
    else
        svarx = varx;
    end
end    
if isnumeric(typ) && ~isempty(typ)
    sthresh = typ;
else
    if strcmpi(typ,'abs') || strcmpi(typ,'dydx') || strcmpi(typ,'polyfit') || strcmpi(typ,'conv')
        styp = typ;
    else
        svarx = typ;
    end
end

if isempty(sthresh)
    sthresh = 0.05;
end
if isempty(svarx)
    svarx = 'Time';
end
if isempty(styp)
    styp = 'dydx';
end

% Identify the right optional variables and defaults
thresh = sthresh;
varx = svarx;
typ = styp;

%% Filter points which are anomalous

% Keep track of the number of points removed
crem = 0;
ctot = 0;

for n=1:length(out)
    
    
    %% Sort the data
    [sx indx] = sort(out(n).data.(varx));
%     out(n).data.(varx) = sx;

    fields = fieldnames(out(n).data);
%     fields(strcmpi(fields,varx) == 1) = [];

    for m = 1:length(fields)
        tt = out(n).data.(char(fields(m)));
        out(n).data.(char(fields(m))) = tt(indx);
    end
    
    
    switch lower(typ)
        case 'abs'
%           perform numerical differentiation          
            dy = diff(out(n).data.(vary));
%             plot(xx,dy)
            ind = abs(dy) > thresh;
            ind = [false; ind];                      % Since diff := x(n) - x(n-1)

            
        case 'dydx'
%           perform numerical differentiation
%           determine mid-points of differential
            xx = out(n).data.(varx)(1:(end-1)) + diff(out(n).data.(varx))/2;
            dydx = diff(out(n).data.(vary))./diff(out(n).data.(varx));
%             plot(xx,dydx)
            ind = abs(dydx) > thresh;
            ind = [false; ind];                      % Since diff := x(n) - x(n-1)
            
        case 'polyfit'
%           fit a second order polynomial to the varx vs vary plot and
%           eliminate points which deviate significantly from the trend
            xx = out(n).data.(varx);
            yy = out(n).data.(vary);
            
            nan_points = [isnan(xx) | isnan(yy)];
            p = polyfit(xx(~nan_points),yy(~nan_points),2);
            yf = polyval(p,xx);
            dif = abs(yy(~nan_points) - yf(~nan_points));
            
%           Determine best suitable threshold to eliminate suitable
%           data points
            [n_i dif_i] = hist(dif(~nan_points),100);
            S_t = 0;
            for m = 1:sum(n_i)
                S_t = S_t + n_i(m)/sum(n_i);
                if S_t > 0.68
                    break
                end
            end

            % Let us consider the deviation of each point away from the
            % intrapolation. If the standard deviation and %68 of data deviate 
            % from the polynomial, it suggests that the errors are random 
            % and can be assumed to be independent. We can then apply only 
            % the user defined threshold on the data
            if std(dif)/dif_i(m) > 0.7 && std(dif)/dif_i(m) < 1.3
                % The threshold is the user-defined threshold
                thresh = sthresh;
            else
                % The threshold is defined by removing ~32% of data
                thresh = dif_i(m);
            end            
            
            % Set the new threshold and determine indicies of bad points
            ind = dif > thresh;
            
            %% Bug testing:
            % -------------------------------------------------------------
            figure, 
            subplot(2,1,1)
            plot(out(n).data.('Temperature2')(ind),out(n).data.('RealPart1a')(ind),'r',...
                 out(n).data.('Temperature2')(~ind),out(n).data.('RealPart1a')(~ind),'g')            
            subplot(2,1,2)
            plot(yy(ind),xx(ind),'r',yy(~ind),xx(~ind),'g')
            % -------------------------------------------------------------
            
        case 'conv'
            % Use 3 point averaging to filter the data
            span = 3; % Size of the averaging window
            window = ones(span,1)/span;
            out(n).data.(vary) = convn(out(n).data.(vary),window,'same');
            
    end
    
    if ~strcmpi(typ,'conv')
        for m = 1:length(fields)
            tt = out(n).data.(char(fields(m)));
            tt(ind) = [];                                  % remove points
            out(n).data.(char(fields(m))) = tt;
        end

        crem = crem + sum(ind);
        ctot = ctot + length(ind);
    else
        crem = 0;
        ctot = 1;
    end
    
end

disp(['Acceptance rate = ' num2str((1 - crem/ctot)*100,'%3.2f') '%'])


end
