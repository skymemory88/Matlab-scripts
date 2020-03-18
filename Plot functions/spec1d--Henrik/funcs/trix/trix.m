function [y, name, pnames, pin]=trix(x,p, flag)
% trix      : 4D convoluted fit function for TAS
% MFIT function [y, name, pnames, pin]=trix(x,p, flag)
% for triple-axis data. A S(Q,w) is convoluted with the
% 4D resolution function, and then a Lorentzian peak is
% added to simulate an incoherent peak.
%
% Calls to: rescal, trixsqw, rc_re2rc, rc_cnmat, rc_savp, mc_monte
%	
% Des McMorrow 25.1.96

global GLOBAL_LENX GLOBAL_PRES GLOBAL_PSCAN GLOBAL_METHOD
global GLOBAL_B_MAT GLOBAL_FWHM GLOBAL_R0

if nargin==2;

%----- Get method from Parameter window
   method = findobj('tag','hrc_rescal_method');
   if ~isempty(method)
	method=get(method,'Userdata');
   else
	return
   end

%----- Rescal params from rescal window 
   pres=rc_savp('respars')';

%----- Scan params from rescal window 
   pscan=rc_savp('tpars')';
   dscan=(pscan(5:8)-pscan(1:4));

%----- Limits on scan
   lx=length(x);
   xstart=x(1);
   xend=x(lx);
   y=zeros(size(x));
   
%----- Read rescal parameters from window

   pres=rc_savp('respars')'; % Rescal:Parameters:sections Spectro + Lattice in row

%----- Read scan parameters from window trixpar window

   pscan=rc_savp('tpars')'; % Rescal:Parameters:section Scan with Mfit in row

   dscan=(pscan(5:8)-pscan(1:4));

%----- Check to see if the scan parameters or spectrometer parameters
%----- have changed - if not use the previously calculated values of
%----- the resolution ellipsoid, otherwise recalculate these quantities. 
compared = 0;
if ~isempty(GLOBAL_LENX)
	compared = compared + sum(lx==GLOBAL_LENX);
end
if ~isempty(GLOBAL_PRES)
	compared = compared + sum(pres==GLOBAL_PRES);
end
if ~isempty(GLOBAL_PSCAN)
	compared = compared + sum(pscan==GLOBAL_PSCAN);
end
if ~isempty(GLOBAL_METHOD)
	compared = compared + sum(method==GLOBAL_METHOD);
end

totlen=(length(lx)+length(pres)+length(pscan)+length(method));

if (compared~=totlen)

%----- Units: f converts from energy units into k^2 (0.4826 for meV, 1.996854 for THz)
%      At the moment, only works for meV!

   f=0.4826;
   mon_flag=p(13);               %  mon_flag: 1=monitor, 0=time

   GLOBAL_LENX=lx;
   GLOBAL_PRES=pres;
   GLOBAL_PSCAN=pscan;
   GLOBAL_METHOD=method;

disp('Rescal : Compute 4D-TAS Resolution function')

%----- Calculate Q2c matrix with Rescal

   Q2c = rc_re2rc( pres(19:21), pres(22:24), pres(31:33) );

%----- Now work out transformations

   A1=pres(25:27)'; % axis 'A' in Lattice in cols
   A2=pres(28:30)'; % axis 'B' 

   V1=Q2c*A1; % in cols
   V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

   V3=cross(V1,V2);
   V2=cross(V3,V1);
   V3=V3/norm(V3);
   V2=V2/norm(V2);
   V1=V1/norm(V1); % in cols

   U=[V1 V2 V3]';

%----- S transformation matrix from (h,k,l) to V1,V2,V3

   S=U*Q2c;     % This is used to bring the CN matrix into a defined frame.

%----- Calculate step size

%   dscan=(pscan(5:8)-pscan(1:4)); % hkle step

   GLOBAL_B_MAT = zeros(lx,16);
   GLOBAL_FWHM = zeros(lx,4);
   GLOBAL_R0 = zeros(lx,1);
                                                
   for j=1:lx

      	pres(31:34)=pscan(1:4)+(x(j)-xstart)/(xend-xstart)*dscan(1:4); % current hkle 

%----- Q vector in cartesian coordinates

      	Qcart=Q2c*pres(31:33)';
      	Qmag=norm(Qcart);
 
%----- Calculate resolution matrix in Q frame. 
%      NP is resolution matrix in Qx, Qy & Qz

      	[R0,NP,vi,vf,Error]=feval(method,f,Qmag,pres,mon_flag);

%----- Use correction of resolution volume
%!!! This part is only a normalisation factor.
%!!! It does not affect the b matrix calculated lower down

      	R0_corrected=R0/(sqrt(det(NP))/(2*pi)^2); % corrected resolution
						% volume as Monte Carlo
						% integral is over a normalised
						% ellipsoid.

%----- Work out angle of Q wrt to V1, V2

      	TT=S*pres(31:33)';
      	cos_theta=TT(1)/norm(TT);
      	sin_theta=TT(2)/norm(TT);

%----- Rotation matrix from Q to V1,V2,V3

      	R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

      	T=zeros(4);
      	T(4,4)=1;
   	T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 

      	M=T'*NP*T; 

	[V,E]=eig(M);
	fwhm=1./sqrt(diag(E))';
	b_mat=reshape(inv((V'))',1,16);

	GLOBAL_B_MAT(j,:)=b_mat; % lx*16
	GLOBAL_FWHM(j,:)=fwhm; % lx*4
	GLOBAL_R0(j)=R0_corrected;

    end % for

    fvcpx = sum(abs(imag(GLOBAL_FWHM')));
    if any(fvcpx)
	disp('Warning : Physical instrumental limits reached.');
	fvrel = find(~fvcpx);
	for j=find(fvcpx)
		[dummy,k]=min(abs(fvrel-j));	% closest element real ellipsoid.
		GLOBAL_FWHM(j,:) = GLOBAL_FWHM(k,:);
		GLOBAL_R0(j) = GLOBAL_R0(k);
		GLOBAL_B_MAT(j,:) = GLOBAL_B_MAT(k,:);
	end
	fprintf(1,'          Fixing %i ellipsoids on scan.\n',length(find(fvcpx)));
	GLOBAL_FWHM = real(GLOBAL_FWHM);
	GLOBAL_R0 = real(GLOBAL_R0);
	GLOBAL_B_MAT = real(GLOBAL_B_MAT);
    end

end % if totlen

%----- Calculate step size  

% pscan : HKLE start, HKLE end
% p : slopes Qx (lin) Qy (quad) Qz (quad), E, bragg point p(5 6 7) Gamma(9) T(11)
	
% get energy at zone center (works for all phonons and Braggs except crossings...). Qx is linear
p(4)=p(4)-p(1)*sqrt((pscan(1)- p(5))^2+(pscan(2)- p(6))^2+(pscan(3) - p(7))^2); %  TA/TO/LA-Phonon  Slope on Qx, zone center gap passed to mex file.  On Qh Qk Ql
%p(4) = p(4)-abs(p(1)*(pscan(2)- p(6)));

% asserting that scan is on phonon polarization vector...
            
NMC=p(12);
for j=1:lx
  pres(31:34)=pscan(1:4)+(x(j)-xstart)/(xend-xstart)*dscan(1:4);
  partemp=[ p(1) p(2) p(3) p(4) p(5) p(6) p(7) p(9) p(11)];
  [smcmex]=mcint(NMC,pres(31:34),...
  	  partemp,GLOBAL_FWHM(j,1:4),GLOBAL_B_MAT(j,1:16));
  y(j)=GLOBAL_R0(j)*smcmex;   
end

%----- Add incoherent peak
y=p(8)*y+p(10)+p(14)*p(15)^2./((x-p(16)).^2+p(15)^2);

else

%----- Mfit parameter window

	y=[];
	name='TRIX';
	mf_msg('Click on Guess button to have some more info...');
	pnames=str2mat('Slope Qx','Slope Qy','Slope Qz','Energy',...
                       'Hzero','Kzero','Lzero',...
                       'Amp','Gamma','Background');
        pnames=str2mat(pnames,'Temperature','NMC','Mon Flag',...
                              'AMPinc','KAPPAinc','CENinc');

	if flag==1, pin=[ 20 0 0 1.05 2 0 0 15000 0.05 1 20 100 1 0 0.001 0 ]; else pin = p; end
	if flag==2
		mf_msg('Sorry, this function does not support guess.')
		disp('TRIX v1 4D-TAS function')
		disp('Dispersion curve (DC) is :');
		disp('linear on scan dir, quadratic on two other base vectors');
	end

%----- TRIX parameter window for resolution and scan parameters

	rescal('trixmode');


end
