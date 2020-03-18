function [data,datastr,scan]=multidata(file,matlabgui,choosedi)
% function [data,datastr,scan]=multidata(file,matlabgui,choosed)
%
% This function enables to load any data
% uses file selector and field chooser
% matlabgui = 1 for GUI window, 0 for none.
% choosed : set of field number to import, can be empty or not precised
%           - for main field (bigger field), * for all

% Author:  EF <manuf@ldv.univ-montp2.fr> 27.06.97
% Description:  General load routine

% uses : looktxt program

data=[];
choosed=[];
scan=[];
datastr=[];
if ~exist('choosedi') choosedi=[]; end
if ~exist('matlabgui') matlabgui=0; end

%===== Load file ============================================

%--------- Open data file ---------------------------------

[fid,m]=fopen(file,'r');
if (fid<0)
   disp([ 'File ' file ' could''nt be opened']);
   disp(m);
   return;
end

filestr = fread (fid);
filestr = setstr(filestr');
fclose(fid);

%------  looktxt analysis ---------------------------------

cmd = [ 'looktxt -F="lktmp001" -t -r -a -f ' file ];
if matlabgui cmd = [ cmd ' -v' ]; end
disp([ 'Starting file analysis :' cmd ]);
eval([ '! ' cmd ]);

patsav = path; % force matlab to rescan directory
addpath(pwd);
path(path);

eval('lktmp001','disp(''No temporary file found.''); return');

path(patsav);
path(path);
delete lktmp001.m;

% now table is 'table'
% [ id type start_pos end_pos rows columns ]
% with type = 1 (num) | 2 (text) | 64 (comment)
if ~exist('table') | isempty(table)
	disp('No data was found in file')
	return
end
test = bitand(1,table(:,2));
numtable = table(find(test),:);
[n,c] = sort(numtable(:,5).*numtable(:,6));
numtable = numtable(c,:);
[n,c] = size(numtable);
numtable = numtable((n:-1:1),:); % revert order : maxid is first
[maxel,maxid] = max(numtable(:,5).*numtable(:,6));
if (n > 30)
	fprintf(1,'Warn : Too many (%i) fields... \n',n);
	fprintf(1,'Press Ctrl-C to abort if needed\n\n');
end

endheaderpos = min([ numtable(maxid,3) numtable(maxid,4) ]);
if (endheaderpos > 500)
	endheaderpos = min([ numtable(1,3) numtable(1,4) ]);
end

header = filestr(max(1,endheaderpos-500):endheaderpos);
fprintf(1,'Header : ...\n%s\n',header);

if (n>1)
  fprintf(1,'Found %i fields\n',n);

  if (matlabgui & isempty(choosedi))


   n=min(n,30);
   fprintf(1,'getting %i biggest fields for GUI\n',n);

  % now create check boxes selector

  %=========== Make checkbox selection window... =======================

  hmf_load=findobj('Tag','mf_cb');
  if ~isempty(hmf_load)
   delete(hmf_load);
  end
   newflag=1;

%------- Create figure window -------------------------------------
   hmf_load=figure('Position',[400 50 174 n*26+50],...
   		'Tag','mf_cb',...
   		'MenuBar','none',...
   		'Name','MFIT: Multi Load',...
   		'Color',get(0,'DefaultUicontrolBackgroundColor'), ...
   		'Resize','off',...
		'Visible','on',...
   		'NumberTitle','off');

%-------- 'Choose field' label ------------------------------------------
   uicontrol(hmf_load,...
   		'Style','text',...
   		'ForegroundColor',[1 1 1],...
   		'Position',[10 (n*26+31) 150 16],...
		'Visible','on',...
   		'String','Choose field');

%--------- check boxes -------------------------------------
   for i=1:n
  	  h(i)=uicontrol(hmf_load,...  
		'Style','checkbox',...       
		'ForegroundColor',[1 1 1],...              
   		'Position',[5 (26*(n+1-i)+2) 150 22],...
		'Visible','on',...
   		'String',sprintf('n%i (%ix%i)',numtable(i,1),numtable(i,5),numtable(i,6)));
   end


%------- ok and cancel buttons -------------------------------------------
   uicontrol(hmf_load,...
   		'Style','PushButton',...
   		'String','Ok',...
   		'Tag','mfcb_ok',...
		'Visible','on',...
   		'Position',[22 4 50 20]);
   uicontrol(hmf_load,...
   		'Style','PushButton',...
   		'String','Cancel',...
   		'Tag','mfcb_cancel',...
   		'CallBack','delete(gcf)', ...
		'Visible','on',...
   		'Position',[82 4 50 20]);

  delete(gca);
  set(hmf_load,'Visible','on');          % Switch window on
  h = findobj('Style','checkbox');
  for i=1:n
	set(h(i),...
	    'String',sprintf('n%i (%ix%i)',numtable(i,1),numtable(i,5),numtable(i,6)));
	if (i == maxid)
		set(h(i),'Value',1);
	end
  end
%============ Now wait until ok or cancel pressed =========================

  drawnow;
  waitforbuttonpress;
  while ~(strcmp(get(gco,'Tag'),'mfcb_ok') | strcmp(get(gco,'Tag'),'mfcb_cancel'))
   	waitforbuttonpress;
	drawnow
  end

%========= Extract relevant data columns ==================================	
  if strcmp(get(gco,'Tag'),'mfcb_ok')
	for i=1:n
	  tmp =  get(h(i),'Value');
	  if (tmp ~= 0)
 	    choosed = [ choosed i ];  
	  end        
	end
  	set(hmf_load,'Visible','off');          % Switch window off
  else
	choosed = [];
	delete(hmf_load);
	return;
  end

  else % no matlabgui or choosed in input args
	for i=1:n
		fprintf(1,'%i - n%i = (%ix%i)\n',i,numtable(i,1),numtable(i,5),numtable(i,6));
	end
	if (isempty(choosedi) | strcmp(choosedi,'-')) choosed=maxid;
	elseif strcmp(choosedi,'*') choosed = 1:n; 
	else choosed = choosedi; end
	if isstr(choosed) choosed = str2num(choosed); end
  end % end if matlabgui
else
	fprintf(1,'Found one field... ok\n');
	choosed = 1;
end

if (isempty(choosed))
	disp('Nothing imported.');
	return;
end

data = [];
n = length(choosed);
fprintf(1,'Getting %i fields : ',n);

nr = sum(numtable(choosed,5));
nc = max(numtable(choosed,6));
data = zeros(nr,nc);
r=1;
for i=1:n
	nr = numtable(choosed(i),5);
	nc = numtable(choosed(i),6);
	ni = ['n' num2str(numtable(choosed(i),1))];
	fprintf(1,'%s  ',ni);
	eval([ 'nf = ' ni ';' ]);
	[nr,nc] = size(nf);
	data((r:(r+nr-1)),1:nc) = nf;
	r=r+nr;
end
[nr,ncol] = size(data);
fprintf(1,'\nData matrix has %ix%i elements.\n',nr,ncol);

datastr='';
choosed=choosed(:)';
scan=num2str(choosed);
