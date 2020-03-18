function [data,datastr]=morpheusheader(file, name)
%
% function [data,datastr]=ntasdata(file)
%
% MATLAB function to load NEW TASCOM file, returning
% the whole  file as an array 'data', with headings 'datastr'
%
% DFM: 27.4.97
%
%----- Open TASCOM data file

data=[]; datastr=[];

dferror=0;

fid=fopen(file,'r');

if (fid <0) 

   datafile=0;
   dferror=1;
   return;
   
end  

%----- Read through data file 

test='zzzz';
data=[];
datastr=[];

while strncmp(test,name,3)==0

    dataline=fgetl(fid);
    if ~isempty(dataline)
	if dataline==-1, dferror=1; return, end 
    end
    test=strtok(dataline);
    if ~isempty(test), test=dataline(1:4); end
           
end

%----- lesen, datastr: name, data: wert

[datastr,dataline]=strtok(dataline); %name
[wort,dataline]=strtok(dataline); %gleichheitszeichen
[wort,dataline]=strtok(dataline); %wert
data=str2num(wort);

fclose(fid);

return
