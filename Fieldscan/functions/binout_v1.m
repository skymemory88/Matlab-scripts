function bout = binout_v1(out,xrange,varx)
% Function to bin the data into bins xrange of variable varx and calculate
% errors based on the standard deviation of the binned data
%
% INPUT:
% out:                  structure containing data
% xrange:               range over which to bin, given as xmin:dx:xmax
% varx:                 parameter which to bin
%
% OUTPUT:
% out:                  binned structure out
%


out = mergeout_v1(out,varx);               % merge all data for binning
%!!! Note that headers, instrument names, filepaths and filenames are not
% merged to retain same structure as input "out"

for n=1:length(out)
    bout(n) = blankstructure;                  % blank structures
end

% Determine the parameters contained in the structure out
fields = fieldnames(out(1).data);

ind = [];

for n = 1:length(fields)
%     if ~isempty(regexpi(fields{n},varx))
%         % Remove the bin parameter from the list
%     else
        if strcmp(fields{n}(1),'d')
            % This variable contains data errors
        else
            % This variable contains data
            ind = [ind; n];
        end
%     end
end

fvar = {fields{ind}}';
flag1 = 0;

% Bin data according to parameter varx grouped into xrange bins
for n = 1:length(out)

    % Select the parameter which determines the bins
    x = out(n).data.(varx);

%     bout.data.(varx) = x;
    
    [h,whichBin] = histc(x, xrange);
    numBins = length(xrange) - 1;
    
    for ii = 1:numBins
        flagBinMembers = (whichBin == ii);
        for m = 1:length(fvar)

            tt = out.data.(char(fvar(m)));
            try
                dtt = out.data.(['d' char(fvar(m))]);
            catch ME
                disp(['oh-oh, missing errors: ' char(fvar(m))])
            end
            
            
            try
                
            bout.data.(char(fvar(m)))(ii,1) = mean(tt(flagBinMembers));
            
                %
                if isempty(dtt)
                    bout.data.(['d' char(fvar(m))])(ii,1) = std(tt(flagBinMembers));                    
                end
                %
                if all(isnan(dtt(flagBinMembers)))
                    % Data errors are missing, use standard deviation of the
                    % parameter as a measure of uncertainty
                    bout.data.(['d' char(fvar(m))])(ii,1) = std(tt(flagBinMembers));
                    flag1 = 1;
                end            
            if all(isnumeric(dtt(flagBinMembers))) && ~all(isnan(dtt(flagBinMembers)))
                % Data errors are missing, use standard deviation of the
                % parameter as a measure of uncertainty
                bout.data.(['d' char(fvar(m))])(ii,1) = norm(dtt(flagBinMembers));
                flag1 = 1;
            end
            catch ME
%                 keyboard
            end
        end
    end


end

end
