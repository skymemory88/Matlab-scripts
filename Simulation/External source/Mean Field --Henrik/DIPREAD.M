function [q,d11,d12]=dipread

fid=fopen('ftn08');

s=fgetl(fid);

q=[];
d11=[];
d12=[];

while ischar(s);

while isempty(findstr('Q',s))
  s=fgetl(fid);
  if ~ischar(s)
    return
  end
end

n=findstr('Q',s);
q=[q;str2num(s(n+4:41))];
s=fgetl(fid);
s=fgetl(fid);
s=fgetl(fid);

s=fgetl(fid);
n1=str2num(s);
s=fgetl(fid);
n2=str2num(s);
s=fgetl(fid);
n3=str2num(s);
d11=[d11;n1(1)+i*n1(4) n2(2)+i*n2(5) n3(3)+i*n3(6)];

s=fgetl(fid);
s=fgetl(fid);

s=fgetl(fid);
n1=str2num(s);
s=fgetl(fid);
n2=str2num(s);
s=fgetl(fid);
n3=str2num(s);
d12=[d12;n1(1)+i*n1(4) n2(2)+i*n2(5) n3(3)+i*n3(6)];


end


