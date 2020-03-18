function [ivst,pvst,fvst,evst]=collate3(iscan,iscanend)

iscanst=iscan;
ivst=[];
pvst=[];
fvst=[];
evst=[];
while iscan <= iscanend
  
   sn=num2str(iscan);
  
   s=loads('specbatch',['upd3.03,X=Theta,Y=Exp_Hutch,M=Monitor,S=' sn] );
   [s,fitdata]=fits(s,'Lorz',[3.0e-4 31.2 0.1 2.0e-4],[1 1 1 1]);
   plot(s);
   pause(1)
   
   st=loads('specbatch',['upd3.03,X=Theta,Y=DegK,S=' sn] );
   [th,temp]=extract(st);
   mtemp=mean(temp);
   
   [ii,eii]=peakt(s);
   pstats=peakm(s,5);
   
   ivst=[ivst; mtemp ii eii];
   pvst=[pvst; mtemp pstats(3)];
   fvst=[fvst; mtemp fitdata.pvals'];
   evst=[evst; mtemp fitdata.evals'];
   
   iscan=iscan+2;
   
 end
