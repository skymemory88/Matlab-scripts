% Program to plot all of the phase diagram data from Feb. 28,
% from which we can construct the phase diagram.

if exist('phase_diag_data.mat') % Should we fit it all?
   load phase_diag_data.mat
else
clear s
s.x_label='Field [T]';
s.y_label='Intensity';
s.datafile='50 mK';
s.yfit=[];

% 50 mK
H50 = [4.35 4.2 4.0 4.1  4.15 4.25];
I50 = [874 876  1163 1002  902 897];
E50 = [ 13   13  15  14 13 13];
s.x=H50;
s.y=I50;
s.e=E50;
s.datafile='50 mK';
s50=spec1d(s);
[s50,f50]=fits(s50,'pow',[1325 4.24 1 874],[1 1 1 1]);

% 65 mK
H65 = [4.0 4.05 4.1 4.1  4.1 4.15 4.2 4.25 4.3 4.35];
I65 = [999 937  870 864  865 778  755 747  733 737];
E65 = [14  13   13  13   9   12   12  12   12  12];
s.x=H65;
s.y=I65;
s.e=E65;
s.datafile='65 mK';
s65=spec1d(s);
[s65,f65]=fits(s65,'pow',[1737 4.21 1.2 740],[1 1 1 1]);


% 100 mK
H100 = [4.0  4.1  4.2  4.3  4.25  4.15];
I100 = [1150 1002 886  869  865   913];
E100 = [15   14   13   13   13    13];
s.x=H100;
s.y=I100;
s.e=E100;
s.datafile='100 mK';
s100=spec1d(s);
[s100,f100]=fits(s100,'pow',[2371 4.22 1.4 868],[1 1 1 1]);


% 195 mK  (average of 190 from CG, 200 from RuO2) 
H195 = [4.0  4.05  4.1  4.15 4.2 4.25  4.3  4.0];
I195 = [1168 1088  998  920  889 873   869  1143];
E195 = [15   14    14   13   13  13    13   15];
s.x=H195;
s.y=I195;
s.e=E195;
s.datafile='195 mK';
s195=spec1d(s);
[s195,f195]=fits(s195,'pow',[2228 4.22 1.34 874],[1 1 1 1]);


% 280 mK
H280 = [4.0  4.1  4.15  4.20 4.25  4.3];
I280 = [1174 976  911   871  888   865];
E280 = [15   14   13    13   13    13 ];
s.x=H280;
s.y=I280;
s.e=E280;
s.datafile='280 mK';
s280=spec1d(s);
[s280,f280]=fits(s280,'pow',[3700 4.2 1.6 874],[1 1 1 1]);


% 370 mK
H370 = [3.9  4.0  4.05  4.1  4.15  4.20 4.25];
I370 = [1254 1129 1054  923  879   867  868];
E370 = [15   15   14    13    13   13    13 ];
s.x=H370;
s.y=I370;
s.e=E370;
s.datafile='370 mK';
s370=spec1d(s);
[s370,f370]=fits(s370,'pow',[949 4.11 0.58 874],[1 1 1 1]);


% 460 mK
H460 = [3.8  3.9  4.0  4.05  4.1  4.15  4.20];
I460 = [1294 1197 1016  903  875   853  856];
E460 = [16   15   14    13    13   13    13 ];
s.x=H460;
s.y=I460;
s.e=E460;
s.datafile='460 mK';
s460=spec1d(s);
[s460,f460]=fits(s460,'pow',[805 4. 0.42 874],[1 1 1 1]);

% 550 mK
H550 = [3.85  3.9  3.95 4.0  4.05  4.1];
I550 = [1141  1029 910  864  841   846];
E550 = [15    14    13    13   13    13 ];
s.x=H550;
s.y=I550;
s.e=E550;
s.datafile='550 mK';
s550=spec1d(s);
[s550,f550]=fits(s550,'pow',[5258 4. 1.7 874],[1 1 1 1]);

% 640 mK This temp is funne, excluded from data.
H640 = [3.80  3.85  3.9  3.95 4.0  3.75  3.75  3.7  3.6];
I640 = [1044  913   863  841  836  957   942   1036 1192];
E640 = [14    13    13    13   13  13    13    14   15];
s.x=H640;
s.y=I640;
s.e=E640;
s.datafile='640 mK';
s640=spec1d(s);
[s640,f640]=fits(s640,'pow',[20 5.2 6 775],[1 1 1 1]);

% 730 mK
H730 = [3.45  3.55  3.65  3.6  3.7  3.75  3.8];
I730 = [1318  1186  912   1032  849  840  841];
E730 = [16    15    13    14    13   13   13 ];
s.x=H730;
s.y=I730;
s.e=E730;
s.datafile='730 mK';
s730=spec1d(s);
[s730,f730]=fits(s730,'pow',[1373 3.66 0.67 874],[1 1 1 1]);

% 820 mK
H820 = [3.4   3.45  3.5   3.55  3.6  3.65];
I820 = [1213  1117  966   863   826  821];
E820 = [15    14    13    13    12   12];
s.x=H820;
s.y=I820;
s.e=E820;
s.datafile='820 mK';
s820=spec1d(s);
[s820,f820]=fits(s820,'pow',[2481 3.57 1.02 874],[1 1 1 1]);

% 910 mK
H910 = [3.25  3.3   3.35  3.4   3.45  3.5   3.55];
I910 = [1295  1238  1119  971   864   840   821];
E910 = [16    15    14    13    13    13    12];
s.x=H910;
s.y=I910;
s.e=E910;
s.datafile='910 mK';
s910=spec1d(s);
[s910,f910]=fits(s910,'pow',[1159 3.4 0.5 874],[1 1 1 1]);

% 1000 mK
H1000 = [3.20  3.25  3.3   3.35  3.4   3.15  3.1];
I1000 = [1148  1022  868   815   819  1126  1198];
E1000 = [15    14    13    12    12    14    15];
s.x=H1000;
s.y=I1000;
s.e=E1000;
s.datafile='1000 mK';
s1000=spec1d(s);
[s1000,f1000]=fits(s1000,'pow',[432 3.25 0.12 874],[1 1 1 1]);

% 1090 mK
H1090 = [2.95  3.05  3.10  3.15  3.2  3.25];
I1090 = [1265  1110  941   852   819   829];
E1090 = [15    14    13    13    12    12];
s.x=H1090;
s.y=I1090;
s.e=E1090;
s.datafile='1090 mK';
s1090=spec1d(s);
[s1090,f1090]=fits(s1090,'pow',[1794 3.16 0.88 874],[1 1 1 1]);

% 1135 mK   --- control 'by hand'
H1135 = [2.7  2.8   2.9  2.95  3.0  3.05  3.10];
I1135 = [1375 1262  1114 941   835   827   832];
E1135 = [16   15    14    13    13    12    12];
s.x=H1135;
s.y=I1135;
s.e=E1135;
s.datafile='1135 mK';
s1135=spec1d(s);
[s1135,f1135]=fits(s1135,'pow',[977 2.96 0.43 824],[1 1 1 1]);

% 1186 mK   --- control 'by hand'
H1186 = [2.7  2.75  2.8   2.85 2.9  2.95];
I1186 = [1226 1132  981   859  856   843];
E1186 = [15   15    14    13    13    13];
s.x=H1186;
s.y=I1186;
s.e=E1186;
s.datafile='1186 mK';
s1186=spec1d(s);
[s1186,f1186]=fits(s1186,'pow',[4239 2.9 1.4 846],[1 1 1 1]);

% 1250 mK   --- control 'by hand'
H1250 = [2.55  2.6  2.65  2.7  2.75  2.8];
I1250 = [1239 1144  1061  882  859   854];
E1250 = [15   15    14    12    13    13];
s.x=H1250;
s.y=I1250;
s.e=E1250;
s.datafile='1250 mK';
s1250=spec1d(s);
[s1250,f1250]=fits(s1250,'pow',[1153 2.7 0.6 874],[1 1 1 1]);


% 1380 mK   --- control 'by hand'
H1380 = [2.35 2.45  2.5   2.55  2.6];
I1380 = [1285 1144  1001  877   879];
E1380 = [15   15    14    13    13];
s.x=H1380;
s.y=I1380;
s.e=E1380;
s.datafile='1380 mK';
s1380=spec1d(s);
[s1380,f1380]=fits(s1380,'pow',[902 2.5 0.4 874],[1 1 1 1]);

% 1420 mK   --- control 'by hand'
H1420 = [2.0  2.10  2.15  2.20  2.25];
I1420 = [1191 1041  923   920   900];
E1420 = [15   14    13    13    13];
s.x=H1420;
s.y=I1420;
s.e=E1420;
s.datafile='1420 mK';
s1420=spec1d(s);
[s1420,f1420]=fits(s1420,'pow',[1083 2.15 0.72 915],[1 1 1 0]);
  save phase_diag_data.mat
end % fitting it all

t=[50 65 100 195 280 370 460 550 730 820 910 1000 1090 1135 1186 1250 1380 1420];
hc=[f50.pvals(2) f50.evals(2)
     f65.pvals(2) f65.evals(2)
     f100.pvals(2) f100.evals(2)
     f195.pvals(2) f195.evals(2)
     f280.pvals(2) f280.evals(2)
     f370.pvals(2) f370.evals(2)
     f460.pvals(2) f460.evals(2)
     f550.pvals(2) f550.evals(2)
     f730.pvals(2) f730.evals(2)
     f820.pvals(2) f820.evals(2)
     f910.pvals(2) f910.evals(2)
     f1000.pvals(2) f1000.evals(2)
     f1090.pvals(2) f1090.evals(2)
     f1135.pvals(2) f1135.evals(2)
     f1186.pvals(2) f1186.evals(2)
     f1250.pvals(2) f1250.evals(2)
     f1380.pvals(2) f1380.evals(2)
     f1420.pvals(2) f1420.evals(2)];
exp=[f50.pvals(3) f50.evals(3)
     f65.pvals(3) f65.evals(3)
     f100.pvals(3) f100.evals(3)
     f195.pvals(3) f195.evals(3)
     f280.pvals(3) f280.evals(3)
     f370.pvals(3) f370.evals(3)
     f460.pvals(3) f460.evals(3)
     f550.pvals(3) f550.evals(3)
     f730.pvals(3) f730.evals(3)
     f820.pvals(3) f820.evals(3)
     f910.pvals(3) f910.evals(3)
     f1000.pvals(3) f1000.evals(3)
     f1090.pvals(3) f1090.evals(3)
     f1135.pvals(3) f1135.evals(3)
     f1186.pvals(3) f1186.evals(3)
     f1250.pvals(3) f1250.evals(3)
     f1380.pvals(3) f1380.evals(3)
     f1420.pvals(3) f1420.evals(3)];
bck=[f50.pvals(4) f50.evals(4)
     f65.pvals(4) f65.evals(4)
     f100.pvals(4) f100.evals(4)
     f195.pvals(4) f195.evals(4)
     f280.pvals(4) f280.evals(4)
     f370.pvals(4) f370.evals(4)
     f460.pvals(4) f460.evals(4)
     f550.pvals(4) f550.evals(4)
     f730.pvals(4) f730.evals(4)
     f820.pvals(4) f820.evals(4)
     f910.pvals(4) f910.evals(4)
     f1000.pvals(4) f1000.evals(4)
     f1090.pvals(4) f1090.evals(4)
     f1135.pvals(4) f1135.evals(4)
     f1186.pvals(4) f1186.evals(4)
     f1250.pvals(4) f1250.evals(4)
     f1380.pvals(4) f1380.evals(4)
     f1420.pvals(4) f1420.evals(4)];

t=[t 1450 1535];
hc=[hc;[1 0];[0 0]];

%%
% Conclusions: critical temperatures and fields.
Tc = [50    65    100   195   280   370   460  550];
Hc = [4.16  4.16  4.17  4.17  4.17  4.125 4.07 3.97];
Tc = [Tc  640   730   820   910   1000  1090  1135  1186];
Hc = [Hc  3.85  3.66  3.55  3.45  3.31  3.17  2.97  2.85];
Tc = [Tc  1250  1380  1420  1430  1535];
Hc = [Hc  2.7   2.55  2.15  1.0   0.0];
%plot(Tc,Hc,'r*-')


% From David Bitko's phase diagram:
Tdb(1) = 0.050000;  Hdb(1) = 50.312;
Tdb(2) = 0.10000 ;  Hdb(2) = 49.268;
Tdb(3) = 0.11400 ;  Hdb(3) = 49.030;
Tdb(4) = 0.15000 ;  Hdb(4) = 47.838;
Tdb(5) = 0.20000 ;  Hdb(5) = 46.324;
Tdb(6) = 0.25000 ;  Hdb(6) = 44.899;
Tdb(7) = 0.30000 ;  Hdb(7) = 43.815;
Tdb(8) = 0.40000 ;  Hdb(8) = 41.991;
Tdb(9) = 0.50000 ;  Hdb(9) = 40.352;
Tdb(10) = 0.60000;  Hdb(10) = 39.359;
Tdb(11) = 0.70000;  Hdb(11) = 37.969;
Tdb(12) = 0.84734;  Hdb(12) = 36.400;
Tdb(13) = 1.0177 ;  Hdb(13) = 33.793;
Tdb(14) = 1.1948 ;  Hdb(14) = 30.240;
Tdb(15) = 1.3321 ;  Hdb(15) = 26.791;
Tdb(16) = 1.4652 ;  Hdb(16) = 22.000;
Tdb(17) = 1.5253 ;  Hdb(17) = 4.9840;
Tdb(18) = 1.5380 ;  Hdb(18) = 0.0000;
Tdb = Tdb.*1000;
Hdb = Hdb./10;

if 1==1
clf

%h=errorbar(t,hc(:,1),hc(:,2),'lino');
h=errorbar(t,hc(:,1),hc(:,2),'o');
set(h,'markersize',5,'linewidth',2')
set(gca,'position',[.15 .15 .8 .8],'fontname','times','fontsize',30)
xlabel('Temperature [mK]')
ylabel('Critical field [T]')
%title('Phase diagram of LiHoF_4, (101) peak intensity')
axis([0 1700 0 5.5])

hold on
  h2=plot(Tdb,Hdb,'r');
hold off
set(h2,'linewidth',2)

legend([h(1) h2],'(1,0,1)','Bitko et al.',3)

print -depsc -loose lihof4_phasediag.eps
end

% Make schematic phase diagram figure (for talks)
if 1==0
   clf
   
   Hc=5.1;
   
   % Shade the quantum critical region.
   t=0:10:1600;
   N=10;
   for n=-N:N
      tt(n+N+1,:)=t;
      hp(n+N+1,:)=n/N*t.^1.2/max(t).^1.2;
      cp(n+N+1,:)=(N-abs(n))/N*fliplr(t);
   end
   pcolor(tt,3*hp+Hc,cp)
   colormap(flipud(hot))
   shading flat
   
   % Plot phase boundary line
   t=[0 Tdb(1:end-2) 1500 Tdb(end-1:end)];
   h=[Hc Hdb(1:end-2) 1.7 Hdb(end-1:end)];
   hold on
     fill([t 0],[h 0],[.2 .6 .8])
   hold off
   h=line(t,h,'color','k','linewidth',2);
   
   
   % And David Bitko's data
   h2=line(Tdb,Hdb,'marker','o','linewidth',2,'color','r','linestyle','none');
   
   nicefig
   set(gca,'fontsize',24)
   xlabel('Temperature [mK]')
   ylabel('Magnetic field [T]')
   
   axis([0 2500 0 8])
   ax=axis;
   ay=ax(3:4)/10;
   ax=ax(1:2)/10;
   
   text(ax*[9.5;.5],ay*[9;1],'Ferromagnet','fontname','times','fontsize',24)
   text(ax*[.5;9.5],ay*[1;9],'Paramagnet','fontname','times','fontsize',24,'verticalal','top','horizontalal','right')
   text(ax*[3;7],ay*[8.5;1.5],str2mat('Thermal','fluctuations'),'fontname','times','fontsize',24)
   text(ax*[9.5;.5],ay*[1;9],str2mat('Quantum','fluctuations'),'fontname','times','fontsize',24)
   text(ax*[6;4],ay*[2.8/8;5.2/8]*10,str2mat('Quantum','critical'),'fontname','times','fontsize',24)
      
   print -depsc -loose lihof4_phasediag_qpt.eps
   print -dtiff -loose lihof4_phasediag_qpt
   
end