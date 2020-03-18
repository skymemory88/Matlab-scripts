t=0.05;
epsilon=0.01;
omega=[0.1:0.05:5];
H=[];
omgs=[];
for h=2:0.2:4
    S=chiqwn([2.01 0 0], omega, epsilon, h, 0.1);
    plot(omega,S)
   % waitforbuttonpress()
    H=[H,h];
    [m,ind]=max(S);
    omgs=[omgs,omega(ind)];
end
plot(H,omgs)