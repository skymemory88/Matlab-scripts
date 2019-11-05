function out = readdata_v3(filepath,filename,run)
%
% =======================================================
% Load data from data file generated by Kelvinox-LQM
%
% INPUT:
% filepath:     directory containing data files
% filename:     filename of the form 'YEAR_MN_'
% run:          file number, ie run = [320 321];
%
% -------------------------------------------------------
%
% OUTPUT:
% out(m):       structure containing data extracted 
%               from the file
%
% =======================================================


fields = {  'Time',...
            'DCField1',...
            'DCField2',...
            'Temperature1',...
            'Temperature2',...
            'Frequency1',...
            'ACField1',...
            'RealPart1a',...
            'ImagPart1a',...
            'RealPart1b',...
            'ImagPart1b',...
            'Frequency2',...
            'ACField2',...
            'RealPart2a',...
            'ImagPart2a',...
            'RealPart2b',...
            'ImagPart2b'};
offset = 0; %compensate is used to make the matrice's dimensions match when not all data are linked properly during the measurement.
nf = length(fields) - offset;
curdir = cd;
tot_time = 0;

for m = 1:length(run)
    mm = [];
    tic
    cd(filepath)

    %file = [filename sprintf('%04.0f',run(m)) '.dat'];
    file = [filename '.dat'];

    fid = fopen(file);
    ll = 1;
    % Add serial time to the data extraction process
    fid = fopen(file);
    tline = importdata(file,',',4);     % The forth line contains string with a time stamp
    date = char(tline(4));  %convert the cell type to string
    t0 = date(7:end);    %select the date-time part of the text
    t0 = datenum(t0, 'dd/mm/yyyy HH:MM:SS');
    fclose(fid);
    
    mm = readtable(file,'Delimiter',{',','\t'});    %import the date from file
    mm = table2array(mm); %convert the data structure to matrix
    
    for n = 1:nf
        out(m).data.(fields{n}) = mm(:,n);   
    end
    
    % Check file for frequency scan from ZVL
    if size(mm,2) > nf     %Counting the elements of each line to check ZVL data
        N = (size(mm,2) - nf)/3;
        out(m).data.ZVLfreq = mm(:,(nf+1):(nf + N));
        out(m).data.ZVLreal = mm(:,(nf+N+1):2:(nf + 3*N));
        out(m).data.ZVLimag = mm(:,(nf+N+2):2:(nf + 3*N));
    end
    
    out(m).headers = '';
    out(m).instrument = '';
    out(m).filepath = filepath;
    out(m).filename = file;
    out(m).run = run(m);
    out(m).data.SerialTime = t0 + out(m).data.Time/86400;
    out(m).starttime = datestr(t0);  % Time when the file was created
    out(m).nop = N;
    
    disp(['...loaded ' out(m).filename ' ' num2str(length(out(m).data.SerialTime)) ' data in ' num2str(toc,3) ' s'])
    tot_time = tot_time + toc;
end

if length(run) > 1
    disp([' Total time: ' num2str(tot_time/60,3) ' min'])
end

cd(curdir)
end

