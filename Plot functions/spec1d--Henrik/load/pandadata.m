function [data,datastr,header,datahead,stepvec,stepstr]=pandadata(filename)
%
% function [data,datastr,header,head,stepvec,stepstr]=pandadata(filename)
%
% MATLAB function to read a PANDA triple axis data file from FRM II
%
% Adapted from illdata by ARW 14.12.05
% Last modified:  ARW 4.01.06
%

%--------- Initialize arrays ------------------------------

data=[];datastr=[];
stepvec=[];stepstr=[];com=[];header=[];
dathead = []; extrasearch = [];
head=[];

%--- temporary 2007/april:
%if exist(filename)~=2 | 1
%    source=['\\Pandasrv\data\2007\',filename(findstr(filename,'prop'):end)];
%    destination='H:\LiHoF4\data';
%    copyfile(source,destination,'f');
%end
%--------- Open data file ---------------------------------

[fid,message]=fopen(filename,'r');
if (fid<0)
   fprintf(1,'ERROR on %s open: %s\n', filename, message);
   return;
end
str='zzzzzz';
while strncmp(str,'scan data:',6)==0
    str=fgetl(fid);
end
%---- Find position for the start of the scan data ---------------------------------

%fpos=ffind(filename,'scan data:');
%if isempty(fpos)
%  %fprintf(1,'%s does not contain any data point.\n', filename);
%  return
%end
%fseek(fid,fpos,'bof');
%r=fgets(fid);
%----- Read in the scan command and create scan information
com=fgetl(fid);
bra=strfind(com,'(');
ket=strfind(com,')');
if ~isempty(bra)
stepstr=deblank(com(1:bra-1));
stepvec=com(bra+1:ket-1);
stepvec(strfind(stepvec,','))=' ';
stepvec=stepvec(~isletter(stepvec));
stepvec=sscanf(stepvec,'%f');
end

%----- Read in the headers. 
%----- These need to be doctored to remove semicolons, hyphens, and to label monitors and detectors
head=fgetl(fid);
%
head1=head;
%
%head=['  ' head(1:findstr(head,';')-1) head(findstr(head,';')+1:end) '  '];
head(head==';')=[];
hyphen=strfind(head,' - ');
if ~isempty(hyphen)
    for i=length(hyphen):-1:1
        head=[head(1:hyphen(i)-1) '-' head(hyphen(i)+3:end)];
    end
end
tagmondet=strfind(head,'tor');
for i=length(tagmondet):-1:1
    if i/2 == floor(i/2)
        j='2'; 
    else 
        j='1';
    end
    head=[head(1:tagmondet(i)+2),j,head(tagmondet(i)+3:end)];
end

%----- Need to change the 'timer' label to 'time'
%timetitle=strfind(head,'time')
%if ~isempty(timetitle)
%    head=[head(1:timetitle+3),head(timetitle+5:end)];
%end;
%new trial with header
head1=[head1(1:findstr(head1,';')-1) head1(findstr(head1,';')+1:end)];
%now second header line
%head2=fgets(fid);
%head2=['  ' head2(1:findstr(head2,';')-1) head2(findstr(head2,';')+1:end) '  '];
head2=head1;

%-----Read data -------------------------------------------
datapos = ftell(fid);
r=fgets(fid);
while length(r) > 2
   semicolon=findstr(r,';');
   if ~isempty(semicolon)
      a=sscanf([r(1:semicolon-1),r(semicolon+1:end)],'%f');
      if ~isempty(a)
      data=[data ; a']; end
   end
   r=fgets(fid);
end

%----Process Header to extract labels ---------------------

datastr = [];
p = 1;
lh = length(head);
while p < lh
  letpos = find(isletter(head(p:lh)));
  if ~isempty(letpos)
    istart = letpos(1)+p-1;
    spapos = find(isspace(head(istart:lh)));
    if ~isempty(spapos)
      iend = spapos(1)+istart-2;
      toadd = head(istart:iend);
      if ~isempty(toadd) & find(isletter(toadd))
        if isempty(datastr)
          datastr = toadd;
        else
          datastr = str2mat(datastr,toadd);
        end
        p = iend+1;
      end
    else
      p = lh;
    end
  else
    p = lh;
  end
end

moretimes=strmatch('time',datastr);
if length(moretimes)>1
    datastr(1,1:4)='timt';
end

%----Process 2. Header to extract units ---------------------

datahead = [];
p = 1;
lh = length(head2);
while p < lh
  letpos = find(isletter(head2(p:lh)));
  if ~isempty(letpos)
    istart = letpos(1)+p-1;
    spapos = find(isspace(head2(istart:lh)));
    if ~isempty(spapos)
      iend = spapos(1)+istart-2;
      toadd = head2(istart:iend);
      if ~isempty(toadd) & find(isletter(toadd))
        if isempty(datahead)
          datahead = toadd;
        else
          datahead = str2mat(datahead,toadd);
        end
        p = iend+1;
      end
    else
      p = lh;
    end
  else
    p = lh;
  end
end

%----- Create header with all the instrument configuration data
fseek(fid,0,-1); % rewind;
%header = fread (fid,fpos-1);
header = setstr(header');


fclose(fid);