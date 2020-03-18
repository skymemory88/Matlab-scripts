function dipolq(qvec,neuerstellen,shells)

%Berechnet diopole() fuer gegebene q werte und speichert diese in
%'dipolwerte'. In Chi_qw wird 'Dpq' geladen. 
if nargin<3
    shells=12;
end
if nargin <2
    neuerstellen=0;
end

N=4; % Number of magnetic atoms in unit cell
Dpqneu=zeros(3,3,N,N,size(qvec,1));

% Calculate dipole sum
for n=1:size(qvec,1)
   Dpqneu(:,:,:,:,n)=dipole(qvec(n,:),shells);
end


if neuerstellen==1
    Dpq=Dpqneu;
    qvektor=qvec;
  %  anz=0;
  %  Dpqneu=zeros(3,3,N,N,size(qvec,1));
  %  qvektorneu=zeros(size(qvec,1),3);
else
    load('dipolwerte', 'Dpq', 'qvektor')
    anz=size(Dpq,5);
    Dtemp=zeros(3,3,N,N,size(qvec,1)+anz);
    Dtemp(:,:,:,:,1:anz)=Dpq; %alte Elemente
    Dtemp(:,:,:,:,(anz+1):(size(qvec,1)+anz))=Dpqneu; %neue Elemente anfügen
    qtemp=zeros(anz+size(qvec,1),3);
    qtemp(1:anz,:)=qvektor;%alte Elemente
    qtemp((1+anz):(anz+size(qvec,1)),:)=qvec;%neue Elemente anfügen
    clear Dpq Dpqneu qvektor
    Dpq=Dtemp;
    qvektor=qtemp;
  %  Dt(:,:,:,:,)=[Dpq(:,:,:,:,)
 %   anz=size(Dpq,5);
 %   Dpqneu=zeros(3,3,N,N,size(qvec,1)+anz);
 %   Dpqneu(:,:,:,:,1:anz)=Dpq;
 %   qvektorneu=zeros(anz+size(qvec,1),3);
 %   qvektorneu(1:anz,:)=qvektor;
 %   clear Dpq qvektor
end

%Dpq=Dpqneu;
%qvektor=qvektorneu;
save('dipolwerte', 'Dpq', 'qvektor')