%%%%%%nullfeld

h=0;
hvec=[0 0 h];
temp=0.05;
momente=[1 0 0
    -1 0 0
    -1 0 0
    1 0 0];
momente=3.4*[0.1 1 0
    -0.1 1 0
    -0.1 -1 0
    0.1 -1 0];
momente=[0.5 0.5 0
    -0.5 0.5 0
    -0.5 0.5 0
    0.5 -0.5 0];
momente=randn(4,3);
omega=[-0.1001:0.001:0.5];
epsilon=0.001;
% with hyperfine
%momente=remf(hvec,temp,momente,0,0.00043421);
%chi0=CHI0_W(hvec,temp,momente,omega,epsilon,0,0.00043421);
% without
 momente=remf(hvec,temp,momente,0,0);
 chi0=CHI0_W(hvec,temp,momente,omega,epsilon,0,0);

q=[0 0 1];
chi=chi_qw(q,chi0);
for nq=1:size(q,1)
    for nw=1:length(omega)
        if temp==0
            Skw(nq,nw)=1/pi* ...
                sum(sum((eye(3)-(q(nq,:)'*q(nq,:))/(q(nq,:)*q(nq,:)')).* ...
                imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
        else
            Skw(nq,nw)= 1/pi/(1-exp(-omega(nw)*11.6/temp))* ...
                sum(sum((eye(3)-(q(nq,:)'*q(nq,:))/(q(nq,:)*q(nq,:)')).* ...
                imag(chi(:,:,nw,nq))));%-chim(:,:,nw,nq))));
        end
    end
end
figure(2)
hold on
plot(omega,Skw,'b')
set(gca,'XLim',[-0.04,0.16],'YLim',[0,2.5e4])


%TEX-file schreiben
fid = fopen('bilder_dispersion\panda1.tex', 'wt');
fprintf(fid,'\\documentclass{article}\n\\usepackage{graphicx}\n\\setlength{\\topmargin}{0pt}\n\\addtolength{\\textwidth}{5cm}\n\\setlength{\\headheight}{0pt}\n\\setlength{\\headsep}{0pt}\n\\addtolength{\\textheight}{5cm}\n\\setlength{\\oddsidemargin}{-1cm}\n\n\\begin{document}\n\n\\section{field along a, H=1 kOe, T= 100 mK}\n\n\\noindent\n');
for n=1:82
    nr=num2str(n);
    fprintf(fid, ['\\includegraphics[width=4.5cm]{afeld', nr,'.eps}\n' ]);
end
fprintf(fid,'\n\\end{document}\n');
fclose(fid);
