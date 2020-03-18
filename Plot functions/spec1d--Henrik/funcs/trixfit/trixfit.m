function [y, name, pnames, pin]=trixfit(x,p,flag)
%
% TRIXFIT function [y, name, pnames, pin]=trixfit(x,p, flag)
% for fitting triple-axis data. 
% 
% S(Q,w) is convoluted with the 4D resolution function.
%
% Des McMorrow, October 2001

%----- Define global variables

global global_trixfit
global first_call
global global_b_mat global_fwhm global_R0

%----- Units: f energy units into k^2 (0.4826 for meV, 1.996854 for THz)
%      At the moment, only works for meV!

if nargin==2

%----- Unpack parameters from global_trixfit

   f=0.4826;

%-----  mon_flag: 1=monitor, 0=time

   mon_flag=global_trixfit.monitor_flag;                                             

%----- Monte Carlo Steps for integration

   NMC=global_trixfit.monte_carlo_samples;

%----- Cross-section

   xsec=global_trixfit.xsec_file;

%----- Background definition

   bkgd_file=global_trixfit.bkgd_file;

%----- Correction file

   corr_file=global_trixfit.corr_file;

%-----  Method 

   method=global_trixfit.resolution_method;

%----- Initialise y

   y=zeros(size(x));

%----- Rescal params

   pres=global_trixfit.pres(:);

%----- Calculate Q2c matrix

   if first_call==1

      [Q2c]=rc_re2rc([pres(19) pres(20) pres(21)], ...
                     [pres(22) pres(23) pres(24)],...
                     [pres(31) pres(32) pres(33)]);

%----- Now work out transformations

      A1=[pres(25) pres(26) pres(27)]';
      A2=[pres(28) pres(29) pres(30)]';

      V1=Q2c*A1;
      V2=Q2c*A2;

%----- Form unit vectors V1, V2, V3 in scattering plane

      V3=cross(V1,V2);
      V2=cross(V3,V1);
      V3=V3/sqrt(sum(V3.*V3));
      V2=V2/sqrt(sum(V2.*V2));
      V1=V1/sqrt(sum(V1.*V1));

      U=[V1';V2';V3'];

%----- S transformation matrix from (h,k,l) to V1,V2,V3

      S=U*Q2c;     % This is used to bring the CN matrix into a defined frame.

%----- Calculate resolution widths for scan etc
 
      for j=1:length(x);

         pres(31)=p(1)*(p(1)<1000)+x(j)*(p(1)>=1000);
         pres(32)=p(2)*(p(2)<1000)+x(j)*(p(2)>=1000);
         pres(33)=p(3)*(p(3)<1000)+x(j)*(p(3)>=1000);
         pres(34)=p(4)*(p(4)<1000)+x(j)*(p(4)>=1000);
        
%----- Calculate focusing curvatures        
        
         if strcmp(method,'rc_popma')
            rho=rc_focus(pres);
            pres(66)=rho(1)*(pres(66)<0);
            pres(67)=rho(2)*(pres(67)<0);
            pres(68)=rho(3)*(pres(68)<0);
            pres(69)=rho(4)*(pres(69)<0);
         end
        
%----- Q vector in cartesian coordinates

         Qcart=Q2c*pres(31:33);
         Qmag=sqrt(sum(Qcart.*Qcart));
 
%----- Resolution matrix in Q frame. NP is resolution matrix in Qx, Qy & Qz

         [R0,NP,vi,vf,Error]=feval(method,f,Qmag,pres,mon_flag);

%----- Use correction of resolution volume

         R0_corrected=R0/(sqrt(det(NP))/((2*pi)^2));

%----- Work out angle of Q wrt to V1, V2

         TT=S*[pres(31) pres(32) pres(33)]';
         cos_theta=TT(1)/sqrt(sum(TT.*TT));
         sin_theta=TT(2)/sqrt(sum(TT.*TT));

%----- Rotation matrix from Q to V1,V2,V3

         R=[cos_theta sin_theta 0; -sin_theta cos_theta 0; 0 0 1];

         T=zeros(4,4);
         T(4,4)=1;
         T(1:3,1:3)=R*S;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 

         M=T'*NP*T; 
         
         [V,E]=eig(M);
         fwhm=zeros(1,4);
         fwhm(1)=1/sqrt(E(1,1));
         fwhm(2)=1/sqrt(E(2,2));
         fwhm(3)=1/sqrt(E(3,3));
         fwhm(4)=1/sqrt(E(4,4));
         b_mat=reshape(inv((V'))',1,16);
         
         global_b_mat(j,:)=b_mat;
         global_fwhm(j,:)=fwhm;
         global_R0(j)=R0_corrected;

      end

   first_call=0;

   end

%----- Calculate function

   for j=1:length(x);

        pres(31)=p(1)*(p(1)<1000)+x(j)*(p(1)>=1000);
        pres(32)=p(2)*(p(2)<1000)+x(j)*(p(2)>=1000);
        pres(33)=p(3)*(p(3)<1000)+x(j)*(p(3)>=1000);
        pres(34)=p(4)*(p(4)<1000)+x(j)*(p(4)>=1000);
        h=pres(31); k=pres(32); l=pres(33); w=pres(34);
        
%----- Evaluate the cross-section and do the Monte Carlo integration

       s=feval(xsec,NMC,[h k l w],p,global_fwhm(j,1:4),global_b_mat(j,1:16));

%----- Calculate normalised intensity

       y(j)=global_R0(j)*s;
                 
   end	

%----- Apply correction factor to calculated intensity, eg for lambda/2 in
%      incident beam

   y=feval(corr_file,x,y);

%----- Add background

   y=y+feval(bkgd_file,x,p);

else
%----- Parameter names

   pnam_file=global_trixfit.pnam_file;
         
   name='trixfit';
   pnames=feval(pnam_file,p);      

end