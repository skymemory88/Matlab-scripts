function[x,y,err,xlab,ylab]=llbload(file)
% Format is
%   qh      qk      ql      en      dh      dk      dl       de    np    m    ki
% ....
%abscisse , comptage  , OHM
% [ X CNTS T ]
 
x=[];y=[];err=[];xlab='';ylab='';
disp(['Loading ' file]);
fid=fopen(file,'r');				
if (fid<0) 
	disp('file not found');
return; 
end;

data=[];		
header='';
text=fgetl(fid);
text=fgetl(fid);
text=fgetl(fid);
text=fgetl(fid);

stext=[0 0];
text=fgetl(fid);
stext=size(text);
while (text>0)
	for i=1:stext(1,2)
		if text(1,i)==','
   		text(1,i)=' ';
		end
	end
data=[data;text];
text=fgetl(fid);
stext=size(text);
end
fclose(fid);
data=str2num(data);
x=data(:,1);
y=data(:,2);
err=sqrt(y);
xlab='Energy';
ylab='Intensity';







