function [data,datastr,pres,pscan,stepvec,stepstr,com,header,head,flip]=illdata(filename)
%
% function [data,datastr,pres,pscan,stepvec,stepstr,header,data_head,flip]=illdata(filename)
%
% MATLAB function to read an ILL triple axis data file
%
% DFM 15.5.96 rev EF 07.07.97 (update)  rev ARW 08.09.98 (polarization)
%

%--------- Initialize arrays ------------------------------

data=[];datastr=[];pres=[];pscan=[];
stepvec=[];stepstr=[];com=[];header=[]; flip = 0;

%--------- Open data file ---------------------------------

fid=fopen(filename,'r');
if (fid<0)
   disp('File not found');
   return;
end

%---- Read column headers ---------------------------------

%fpos=ffind(filename,'DATA_:');
% THis line was put in may 99 to read data from v2-flex at HMI.
fpos=ffind(filename,'PNT')-1;
fseek(fid,fpos,'bof');
r=fgets(fid);
head=fgets(fid);
head=['  ' head '  '];

%-----Read data -------------------------------------------
datapos = ftell(fid);
r=fgets(fid);
while length(r) > 2
   r=strrep(r,'*',' ');
   a=sscanf(r,'%f');
   data=[data ; a'];
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

fpos=ffind(filename,'COMND:');
fseek(fid,fpos,'bof');
com=fgets(fid);

fseek(fid,0,-1); % rewind;
header = fread (fid,datapos-1);
header = setstr(header');


%----- Now read the resolution parameters

if nargout > 3
                
   pres=1:42;

   fpos=findstr(header,'PARAM');
   if isempty(fpos)
%	disp('No param section.');
   else
	section=header(fpos(1):length(header));	
   	tosearch=str2mat('DM=','DA=',...
		'ETAM=','ETAA=','ETAS=',...
		'SM=','SS=','SA=',...
		'KFIX=','FX=',...
		'ALF1=','ALF2=','ALF3=','ALF4=',...
		'BET1=','BET2=','BET3=','BET4=',...
		'AS=','BS=','CS=',...
		'AA=','BB=','CC=',...
		'AX=','AY=','AZ=',...
		'BX=','BY=','BZ=');
	[n,c]=size(tosearch);
	x=1:n;
   	for j=1:n
		item=deblank(tosearch(j,:));
		fpos=findstr(section,item);
		if isempty(fpos)
			x(j) = NaN;
		else
			fpos=fpos(1)+length(item);
			x(j) = sscanf(section(fpos:length(section)),'%f');
		end
	end
	pres(1:30)=x;
   end

   fpos=findstr(header,'POSQE');
   if isempty(fpos)
%	disp('No posqe section.');
   else
	section=header(fpos(1):length(header));	
   	tosearch=str2mat('QH=','QK=','QL=','EN=');
	[n,c]=size(tosearch);
	x=1:n;
   	for j=1:n
		item=deblank(tosearch(j,:));
		fpos=findstr(section,item);
		if isempty(fpos)
			x(j) = NaN;
		else
			fpos=fpos(1)+length(item);
			x(j) = sscanf(section(fpos:length(section)),'%f');
		end
	end
	pres(31:34)=x;
   end

   x=[];
   fpos=ffind(filename,'STEPS');
   fseek(fid,fpos,'bof');
   r=fgets(fid);
   steps = [ ' ' r ' ' ];
   lsteps = length(steps);
   stepstr2 = '';

   dqh = 0;
   id = findstr(steps,'DQH');
   if (~isempty(id))
      id1=findstr(steps((id(1)+3):lsteps),'=');
      x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
      dqh = x(1);
      if isempty(stepstr2)
	stepstr2 = 'DQH';
      else
	stepstr2 = str2mat(stepstr2,'DQH');
      end
	stepvec = [ stepvec ; x(1) ];
   end

   dqk = 0;
   id = findstr(steps,'DQK');
   if (~isempty(id))
      id1=findstr(steps((id(1)+3):lsteps),'=');
      x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
      dqk = x(1);
      if isempty(stepstr2)
	stepstr2 = 'DQK';
      else
	stepstr2 = str2mat(stepstr2,'DQK');
      end
	stepvec = [ stepvec ; x(1) ];
   end

   dql = 0;
   id = findstr(steps,'DQL');
   if (~isempty(id))
      id1=findstr(steps((id(1)+3):lsteps),'=');
      x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
      dql = x(1);
      if isempty(stepstr2)
	stepstr2 = 'DQL';
      else
	stepstr2 = str2mat(stepstr2,'DQL');
      end
	stepvec = [ stepvec ; x(1) ];
   end

   den = 0;
   id = findstr(steps,'DEN');
   if (~isempty(id))
      id1=findstr(steps((id(1)+3):lsteps),'=');
      x = sscanf(steps((id(1)+id1(1)+3):lsteps),'%f');
      den = x(1);
      if isempty(stepstr2)
	stepstr2 = 'DEN';
      else
	stepstr2 = str2mat(stepstr2,'DEN');
      end
	stepvec = [ stepvec ; x(1) ];
   end


   pres(35:38)=[ dqh dqk dql den];
   pres(39:42)=[ 0 0 1 1];

   if isempty(stepvec)
   streq=findstr('=',r);
   for i=1:length(streq)
      stepvec =[stepvec sscanf(r(streq(i)+1:length(r)),'%f',1)];
   end
   end

%----Process STEPS header to extract labels ---------------------
   if ~isempty(stepstr2)
	stepstr = stepstr2;
   else
	p=findstr(steps,' ');
	K=p(filter([1 -1],1,p)>1); % start of non blank string analysing differences
	stepstr=[];
	for i=1:3:length(K)-1
		addstr = findstr(steps((K(i)+2):length(steps)),' ');
		addstr = min(addstr(1)+(K(i)+2),lsteps);
		addstr=steps((K(i)+1):addstr);
		if (i==1) stepstr=addstr;
		else
			stepstr=str2mat(stepstr,addstr);
		end
	end
   end
end

%----- Find if polarisation was used

ipol=findstr(head,'PAL');
if ~isempty(ipol)
	flip=1;
else
	flip=0;
end

%----- Calculate the scan limits

if nargout > 4

   pscan(1:4)=pres(31:34);
   [npts,dummy]=size(data);
   pscan(5:8)=pres(31:34)+pres(35:38)*(npts-1);		

end

fclose(fid);
