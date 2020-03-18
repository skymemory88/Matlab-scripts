function Skw=scatt(q,omega,epsilon,t,h)

chi=chi_qw(q,omega,epsilon,t,h);
chim=chi_qw(q,-omega,epsilon,t,h);

Skw=zeros(length(omega),size(q,1));
for nw=1:length(omega)
for nq=1:size(q,1)
  if t>0
    Skw(nw,nq)=1/pi/(1-exp(-omega(nw)*11.6/t))* ...
      sum(sum((eye(3)-(q(nq,:)'*q(nq,:))/(q(nq,:)*q(nq,:)')).* ...
              imag(chi(:,:,nw,nq)-chim(:,:,nw,nq))));
  else
    Skw(nw,nq)=1/pi* ...
      sum(sum((eye(3)-(q(nq,:)'*q(nq,:))/(q(nq,:)*q(nq,:)')).* ...
              imag(chi(:,:,nw,nq)-chim(:,:,nw,nq))));
  end
end
end