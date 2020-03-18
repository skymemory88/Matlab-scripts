function [JEr,altJEr,JErhyp,altJErhyp,JHo,altJHo,history]=LiErxHoyY1_x_yF4(x,y,temp,field,exEr,exHo,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho)
%Compute the magnetization and the alternated magnetization for a range of temperatures and
%fields. (1-x) is the dilution. temp and h are the temperature and field
%arrays. xHypIso is the proportion of nuclear moments isotops for Er. alpha
%is, in the case the sample is an ellipsiod (a,c,c), the ratio a/c (for
%demagnetization tensor.

global strategies;

t=tic; %start a stopwatch timer

mom_Er=[1 0 0
        -1 0 0
        -1 0 0
        1 0 0];

mom_Ho=[0 0 1
        0 0 1
        0 0 1
        0 0 1]*2.6;%*0.0001;%*0.1365;

alt=[1 1 1
    -1 -1 1
    -1 -1 1
     1 1 1];

altJEr=zeros([3,length(temp),size(field,2)]);
altJErhyp=zeros([3,length(temp),size(field,2)]);
altJHo=zeros([3,length(temp),size(field,2)]);
JEr=zeros([3,length(temp),size(field,2)]);
JErhyp=zeros([3,length(temp),size(field,2)]);
JHo=zeros([3,length(temp),size(field,2)]);

history=cell([length(temp),size(field,2)]);

for i = 1:size(field,2)
    for j = 1:length(temp)
        [mom_Er,mom_Er_hyp,mom_Ho,evolution]=remf(field(:,i)',temp(j),mom_Er,mom_Ho,exHo,exEr,x,y,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho);
        altJEr(:,j,i)=mean(alt.*mom_Er);
        altJErhyp(:,j,i)=mean(alt.*mom_Er_hyp);
        JEr(:,j,i)=mean(mom_Er);
        JErhyp(:,j,i)=mean(mom_Er_hyp);
        JHo(:,j,i)=mean(mom_Ho);
        altJHo(:,j,i)=mean(alt.*mom_Ho);
        history{j,i}=evolution;
        if strategies.powerlaw && j>3 && j<length(temp) % Guess starting values for next temperature using powerlaw
            beta=0.5; % assume exponent 1/2. Could generalize to use more points and flexible exponent
            ibeta=1/beta;
            JHo1=squeeze(history{j,i}.ho(end,:,:));
            JHo2=squeeze(history{j-1,i}.ho(end,:,:));
            JEr1=squeeze(history{j,i}.er(end,:,:));
            JEr2=squeeze(history{j-1,i}.er(end,:,:));
            TcHo=((JHo1./JHo2).^ibeta*temp(j-1)-temp(j))./((JHo1./JHo2).^ibeta-1); % can generalize to use more points
            TcEr=((JEr1./JEr2).^ibeta*temp(j-1)-temp(j))./((JEr1./JEr2).^ibeta-1); % can generalize to use more points
            if isempty(find(((TcHo<=0) + (TcHo > temp(end))),1)) % only do this if Tc guess is reasonable
                Aa=JHo2.^ibeta./(TcHo-temp(j-1));
                if isempty(find(~isfinite(Aa)+isnan(Aa),1)) % problem if Tc happens to be equal temp(j-1)
                    JHonew=real((Aa.*(TcHo-temp(j+1))).^beta).*(temp(j+1)<TcHo);
                    mom_Ho=JHonew;
                end
            end
            if isempty(find(((TcEr<=0) + (TcEr > temp(end))),1)) % only do this if Tc guess is reasonable
                Aa=JEr2.^ibeta./(TcEr-temp(j-1));
                if isempty(find(~isfinite(Aa)+isnan(Aa),1)) % problem if Tc happens to be equal temp(j-1)
                    JErnew=real((Aa.*(TcEr-temp(j+1))).^beta).*(temp(j+1)<TcEr);
                    mom_Er=JErnew;
                end
            end
        end
    end
end

altJEr=squeeze(altJEr);
altJErhyp=squeeze(altJErhyp);
altJHo=squeeze(altJHo);
JEr=squeeze(JEr);
JErhyp=squeeze(JErhyp);
JHo=squeeze(JHo);
toc(t)
    
