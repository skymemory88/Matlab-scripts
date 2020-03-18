function chiqw_erstellen()
%clear all
qv=[1:0.025:2]'*[1 0 0];
ow=[0:0.01:0.8];

bild=zeros(3,3,length(ow),size(qv,1));
y=zeros(4,size(qv,1));
%tw=zeros(length(omega),1);
for hn=1:4
    bild=chi_qw(qv,ow,0.01,2+hn,0.4);
    for m=1:size(qv,1)
        for n=1:length(ow)
            tw(n)=norm(imag(bild(:,:,n,m)));
        end
        [wert,indx] = max(tw); 
        y(hn,m)=ow(indx);
    end
end

plot(qv(:,1),y(1,:), qv(:,1), y(2,:), qv(:,1) ,y(3,:),qv(:,1) ,y(4,:))
