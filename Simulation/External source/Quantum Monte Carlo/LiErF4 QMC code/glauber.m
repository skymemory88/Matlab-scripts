function acc=glauber(dE, temp)

if(temp==0)
    if(dE<0)
        acc=1;
    else acc=0;
    end
else
    proba=1/(1+exp(dE*11.6/temp));
    alpha=rand(1);
    if(alpha<=proba)
        acc=1;
    else acc=0;
    end
end