load pg1 pt ph

epsilon=0.001;
omega=[-0.0101:0.0002:0.02];
H=[];
omgs=[];
for qh=1.61:0.05:2.01
    S=chiqwn([qh 0 0], omega, epsilon, ph(2), pt(2));
    plot(omega,S)
    pause(0.1)
    hold on
end
