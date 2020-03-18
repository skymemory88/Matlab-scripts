function [data,datastr,head,headstr]=hmidata(file)
%
% function [data,datastr]=hmidata(file)
%
% MATLAB function to load HMI carres data file, returning
% the whole  file as an array 'data', with headings 'datastr'
%
% HMR 16.8.98

% Initialize output variables
data=[]; datastr=[];

% Open file
dferror=0;
fid=fopen(file,'r');
if (fid <0) 
   datafile=0;
   dferror=1;
   return;
end  

%----- Read through data file 
dataline=getline(fid)
[column_lab,dataline]=strtok(dataline);
while ~isempty(dataline)
   [column_lab,dataline]=strtok(dataline);
   if num2str(column_lab(1))==0
       if column_lab(1)==0
           data=[data, 0,0]
       else
           
   datastr=strvcat(datastr,column_lab);
end


return
dataline='zzzz';
while strncmp(dataline,' NR.',4)==0
  dataline=fgetl(fid);
  if ~isempty(dataline)
	 if dataline==-1, dferror=1; return, end 
  end
end
%----- Create column headers

%[kill,dataline]=strtok(dataline);
while ~isempty(dataline)
   [column_lab,dataline]=strtok(dataline);
   datastr=strvcat(datastr,column_lab);
end
%----- leerzeile eberspringen (Nur bei 2D-Detektor,sonst auskommentieren)
fgetl(fid);
%----- Read data

dataline=str2num(fgetl(fid));
while ~isempty(dataline)
  dataline=[dataline zeros(1,size(data,2)-length(dataline))];
  data=[data;dataline];
  dataline=str2num(fgetl(fid));
end

fclose(fid);

return
