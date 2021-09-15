function sample = LiErxHoyY1_x_yF4_mod(sample)
%Compute the magnetization and the alternated magnetization for a range of temperatures and
%fields. (1-x) is the dilution. temp and h are the temperature and field
%arrays. xHypIso is the proportion of nuclear moments isotops for Er. alpha
%is, in the case the sample is an ellipsiod (a,c,c), the ratio a/c (for
%demagnetization tensor.

% [JEr,altJEr,JErhyp,altJErhyp,JHo,altJHo,Jmom,history]=...
% LiErxHoyY1_x_yF4_mod(x,y,temp,field,exEr,exHo,ErHyp,HoHyp,demagn,alpha,renorm_Er,renorm_Ho)

temp = sample.temp;
field = sample.fields;

global strategies;

t=tic;
sample.mom_Er=[1 0 0
        -1 0 0
        -1 0 0
        1 0 0];

sample.mom_Ho=[0 0 1
    0 0 1
    0 0 1
    0 0 1]*2.6;%*0.0001;%*0.1365;

alt=[1 1 1
    -1 -1 1
    -1 -1 1
    1 1 1];

sample.altJEr=zeros([3,length(temp),size(field,2)]);
sample.altJErhyp=zeros([3,length(temp),size(field,2)]);
sample.altJHo=zeros([3,length(temp),size(field,2)]);
sample.JEr=zeros([3,length(temp),size(field,2)]);
sample.JErhyp=zeros([3,length(temp),size(field,2)]);
sample.JHo=zeros([3,length(temp),size(field,2)]);
sample.Jmom=zeros(4,3,length(temp),size(field,2));

sample.Ham_hyp_Er_TH = zeros(4,length(temp),size(field,2),(2*sample.J_Er + 1)*(2*sample.I_Er + 1),(2*sample.J_Er + 1)*(2*sample.I_Er + 1));
sample.Ham_Er_TH = zeros(4,length(temp),size(field,2),2*sample.J_Er + 1,2*sample.J_Er + 1);
sample.Ham_hyp_Ho_TH = zeros(4,length(temp),size(field,2),(2*sample.J_Ho + 1)*(2*sample.I_Ho + 1),(2*sample.J_Ho + 1)*(2*sample.I_Ho + 1));
sample.Ham_Ho_TH = zeros(4,length(temp),size(field,2),2*sample.J_Ho + 1,2*sample.J_Ho + 1);

sample.history=cell([length(temp),size(field,2)]);

for i = 1:size(field,2)
    for j = 1:length(temp)
        [mom_Er,mom_Er_hyp,mom_Ho,mom_mean,sample,evolution]=...
            remf_mod(sample,field(:,i)',temp(j));
        sample.altJEr(:,j,i)=mean(alt.*mom_Er);
        sample.altJErhyp(:,j,i)=mean(alt.*mom_Er_hyp);
        sample.JEr(:,j,i)=mean(mom_Er);
        sample.JErhyp(:,j,i)=mean(mom_Er_hyp);
        sample.JHo(:,j,i)=mean(mom_Ho);
        sample.altJHo(:,j,i)=mean(alt.*mom_Ho);
        sample.Jmom(:,:,j,i) = mom_mean;
        
        sample.Ham_hyp_Er_TH(:,j,i,:,:) = sample.Ham_hyp_Er;
        sample.Ham_Er_TH(:,j,i,:,:) = sample.Ham_Er;
        sample.Ham_hyp_Ho_TH(:,j,i,:,:) = sample.Ham_hyp_Ho;
        sample.Ham_Ho_TH(:,j,i,:,:) = sample.Ham_Ho;
        
        sample.history{j,i}=evolution;
        if strategies.powerlaw && j>3 && j<length(temp) % Guess starting values for next temperature using powerlaw
            beta=0.5; % assume exponent 1/2. Could generalize to use more points and flexible exponent
            ibeta=1/beta;
            JHo1=squeeze(sample.history{j,i}.ho(end,:,:));
            JHo2=squeeze(sample.history{j-1,i}.ho(end,:,:));
            JEr1=squeeze(sample.history{j,i}.er(end,:,:));
            JEr2=squeeze(sample.history{j-1,i}.er(end,:,:));
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

% altJEr=squeeze(altJEr);
% altJErhyp=squeeze(altJErhyp);
% altJHo=squeeze(altJHo);
% JEr=squeeze(JEr);
% JErhyp=squeeze(JErhyp);
% JHo=squeeze(JHo);
toc(t)
    
