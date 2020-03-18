plot(Skw(5,:))
mwechsel=[];
for nw=1:length(ow)
    maufab=0;
    for nh=2:length(hvek)
        if maufab==1
            if Skw(nw,nh)<=Skw(nw,nh-1)
                maufab=0;
                mwechsel=[mwechsel;[hvek(nh) ow(nw)]];
            end
        else
            if Skw(nw,nh)>=Skw(nw,nh-1)
                maufab=1;
            end
        end
    end
end
mwechsel;
figure 
plot(mwechsel(:,1),mwechsel(:,2))