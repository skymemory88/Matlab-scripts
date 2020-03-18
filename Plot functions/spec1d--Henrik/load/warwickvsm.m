function [data,datastr,parval,parstr,title]=warwickvsm(file)
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

dataline=fgetl(fid); % Skip the line with the date.
disp(['Reading ',file,' from ',dataline])
dataline=fgetl(fid); % Get the parameters for this sweep.
[title,rest]=strtok(dataline,',');
parstr=[];
parval=[];
while ~isempty(rest)
  [partxt,rest]=strtok(rest,',');
  [nparstr,nparval]=strtok(partxt,'=');
  nparstr=deblank(nparstr);
  [a,b]=strtok(nparval);
  [a,b]=strtok(b);
  nparval=str2num(a);
  nparstr=[nparstr,' ',deblank(b)];
  parstr=str2mat(parstr,nparstr);
  parval=[parval;nparval];
end
parstr(1,:)=[];

dataline1=fgetl(fid); % Get the data names
dataline2=fgetl(fid); % Get the data units
[a,rest1]=strtok(dataline1,',');
[b,rest2]=strtok(dataline2,',');
a=fliplr(deblank(fliplr(a)));
b=fliplr(deblank(fliplr(b)));
datastr=[deblank(a),' [',deblank(b),']'];
while ~isempty(rest1) & ~isempty(rest2)
  [a,rest1]=strtok(rest1,',');
  [b,rest2]=strtok(rest2,',');
  a=fliplr(deblank(fliplr(a)));
  b=fliplr(deblank(fliplr(b)));
  datastr=str2mat(datastr,[deblank(a),' [',deblank(b),']']);
end

dataline=fgetl(fid); % Skip numerical format information

data=[];
dataline=fgetl(fid); % Get the data
while dataline~=-1
  data=[data;str2num(dataline)];
  dataline=fgetl(fid);
end

fclose(fid);

return

%----- Create column headers

[kill,dataline]=strtok(dataline);
while ~isempty(dataline)

   [column_lab,dataline]=strtok(dataline);
   datastr=strvcat(datastr,column_lab);

end

%----- Read data

while(1)

   dataline=fgetl(fid);
   if dataline==-1, 
     fclose(fid)
     return 
   elseif ~isempty(dataline)
      a=sscanf(dataline,'%f'); data=[data; a'];
   end

end

fclose(fid)

return
