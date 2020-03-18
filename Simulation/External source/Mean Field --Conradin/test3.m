hvek=[0:0.05:6];
ow=[0.01:0.01:0.5]; %omega
Hcf=cf();
chi0=zeros(3,3,length(ow),length(hvek));
t=0.8;

for jh=1:length(hvek)
    [e,v]=mfieldz(hvek(jh),t,1,0,1,Hcf,0);
    for jw=1:length(ow)
        chi0(:,:,jw,jh)=chi0_w([hvek(jh) 0 0],t,ow(jw),0.01,e,v); 
    end
end
figure
pcolor(hvek,ow,imag(squeeze(chi0(2,2,:,:))))
xlabel('H_x / Tesla')
ylabel('\omega / meV')