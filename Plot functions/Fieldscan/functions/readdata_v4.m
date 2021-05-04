
function out = readdata_v4(file,nZVL)
%
% =======================================================
% Load data from data file generated by the new labview 
% measurement suite
%
% INPUT:
% filepath:     directory containing data files
% filename:     filename of the form 'YEAR_MN_####'
%
% -------------------------------------------------------
%
% OUTPUT:
% out:       structure containing data extracted 
%               from the file
%
% =======================================================


fields = {  'Time',...
            'DCField1',...
            'DCField2',...
            'Temperature1',...
            'Temperature2',...
            'Power1',...
            'Attenuation1',...
            'Custom',...
            'lockin1',...
            'lockin2'
          };
offset = 0; %compensate is used to make the matrice's dimensions match when not all data are linked properly during the measurement.
nf = int16(length(fields) - offset);
curdir = cd;

tic

fid = fopen(file);
tline = importdata(file,',',4);     % The forth line contains string with a time stamp
date = char(tline(4));  %convert the cell type to string
t0 = date(7:end);    %select the date-time part of the text
t0 = datenum(t0, 'dd/mm/yyyy HH:MM:SS');
fclose(fid);

mm = readtable(file,'Delimiter',{',','\t'});    %import the date from file
mm = table2array(mm); %convert the data structure to matrix

for n = 1:nf
    out.data.(fields{n}) = mm(:,n);   
end

% Check file for frequency scan from ZVL
if size(mm,2) > nf    %Counting the elements of each line to check ZVL data
        N = int16(size(mm,2) - nf - 1)/(1+(2*nZVL));
%         N = (size(mm,2) - nf)/(2*(1+nZVL));
        out.data.ZVLfreq = mm(:,(nf+1):(nf + N));
        out.data.ZVLreal = mm(:,(nf+N+1):2:(nf + 3*N)); % S11 data
        out.data.ZVLimag = mm(:,(nf+N+2):2:(nf + 3*N));
    if nZVL == 2
%         out.data.ZVLreal2 = mm(:,(nf+3*N+1):2:(nf + 5*N)); % S21 data
%         out.data.ZVLimag2 = mm(:,(nf+3*N+2):2:(nf + 5*N));
        out.data.ZVLreal2 = mm(:,(nf+4*N+1):2:(nf + 6*N)); % S21 data
        out.data.ZVLimag2 = mm(:,(nf+4*N+2):2:(nf + 6*N));
    end
end

out.headers = '';
out.instrument = '';
[out.filepath, out.filename] = fileparts(file);
out.data.SerialTime = t0 + out.data.Time/86400;
out.starttime = datestr(t0);  % Time when the file was created
out.nop = N;

disp(['...loaded ' out.filename ' ' num2str(length(out.data.SerialTime)) ' data in ' num2str(toc,3) ' s'])
cd(curdir)
end

