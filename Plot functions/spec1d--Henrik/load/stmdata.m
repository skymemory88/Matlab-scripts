function [data,datastr,par]=stmdata(file)
%
% function [data,datastr]=stmdata(file)
%
% MATLAB function to load STM spectra saved from STM-user written by Chris Renner
%
% HMR: 13th of November 2002
%

dferror=0;

fid=fopen(file,'r');
if (fid <0) 
   datafile=0;
   dferror=1;
   return;   
end  

%----- Read through data file 

dataline='zzzz';
data=[];
datastr=[];

dataline=fgetl(fid);
dataline=fgetl(fid);
dataline=fgetl(fid);
dataline=fgetl(fid);
dataline=[fgetl(fid) ' ' fgetl(fid)];
par(1,6)=datenum(dataline)*24*3600;
while strncmp(dataline,'Current nA:',11)==0
    dataline=fgetl(fid);
end
par(1,3)=str2num(fgetl(fid));
while strncmp(dataline,'Bias V:',7)==0
    dataline=fgetl(fid);
end
par(1,2)=str2num(fgetl(fid));
while strncmp(dataline,'xscanner',8)==0
    dataline=fgetl(fid);
end
par(1,4:5)=str2num(fgetl(fid));
while strncmp(dataline,'T:',2)==0
    dataline=fgetl(fid);
end
tmp=str2num(strtok(dataline(3:end)));
if ~isempty(tmp)
   par(1,1)=tmp;
else
   par(1,1)=0;
end

%----- Create column headers

while strncmp(dataline,'DATA::',6)==0
    dataline=fgetl(fid);
end
datastr=fgetl(fid);
datastr=strvcat(datastr,fgetl(fid));
tmp=fgetl(fid);
%----- Read data
while(1)
   dataline=fgetl(fid);
   if dataline==-1, 
      fclose(fid);
      return 
   elseif ~isempty(dataline)
      a=sscanf(dataline,'%f'); data=[data; a'];
   end
end

fclose(fid);
return
