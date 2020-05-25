function measurement(jobid,nummeas,type,eq)
% Function to compute field and temperature scans based on output in files:
% ...\MC\data\LiErF4\tempscans
% INPUT:
% jobid             # of configuration
% nummeas           # of measurement on this configuration
% type:             temp or field scan
% eq:               0 if you already computed the equilibrium and 1 if want
%                   to compute it

if(strcmp(type,'temp'))
    if(eq==1)
        MCscriptEQ(jobid)
    end
    MCmeasuretempscan(jobid,nummeas)
elseif(strcmp(type,'field'))
    if(eq==1)
        MCscriptEQ(jobid)
    end
    MCmeasurefieldscan(jobid,nummeas)   
    
else
    fprintf(stderr,'type must be temp or field !');
    return
end
end