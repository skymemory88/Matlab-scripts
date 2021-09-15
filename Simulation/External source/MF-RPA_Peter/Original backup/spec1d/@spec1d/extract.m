function [varargout]=extract(s1)
%
% function [{x,y,e,yfit}]=extract(s1)
%
% @SPEC1D/EXTRACT Extracts raw data from spec1d object s1, or from
% spec1d array if all objects in array have same number of data points
%
% Version 2.0, February 2001
% Des McMorrow and Henrik Ronnow

nout=nargout;

if ~isa(s1,'spec1d')

   disp('Extract error: input must be a spec1d data object')
   return

end

iobj=length(s1);
lenobj=length(s1(1).x);
x=zeros(lenobj,iobj); 
y=zeros(lenobj,iobj); 
e=zeros(lenobj,iobj); 
yfit=zeros(lenobj,iobj); 

for il=1:iobj

   x(:,il)=s1(il).x;
   y(:,il)=s1(il).y;
   e(:,il)=s1(il).e;
   if isempty(s1(il).yfit)
      yfit(:,il)=-1*ones(size(s1(il).e));
   else
      yfit(:,il)=s1(il).yfit;      
   end   
      
end

if nout>=1
   varargout(1)={x};
end   
if nout>=2
   varargout(2)={y};
end
if nout>=3
   varargout(3)={e};
end
if nout==4
   varargout(4)={yfit};
end