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

dataline='zzzz';
while strncmp(dataline,'READ',4)==0
  dataline=fgetl(fid);
  if ~isempty(dataline)
	 if dataline==-1, dferror=1; return, end 
  end
end
dataline2=fgetl(fid);
if dataline2(1)==' ' | dataline2(1)=='-'
   dataline=[dataline(1:end-2) dataline2(1:end)];
end
% Extract header information
dataline=dataline(5:end);
headstr=[];
head=[];
while ~isempty(dataline)
[t,dataline]=strtok(dataline,'= ');
headstr=str2mat(headstr,t);
[t,dataline]=strtok(dataline,'= ');
head=[head;str2num(t)];
end
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

%----- Read data

dataline=str2num(fgetl(fid));
while ~isempty(dataline)
  dataline=[dataline zeros(1,size(data,2)-length(dataline))];
  data=[data;dataline];
  dataline=str2num(fgetl(fid));
end

fclose(fid);

return
