c     Fortran program to calculate the
c     Fourrier transforms of the dipole-dipole interaction
c     using Ewald's methos as refined by
c     Bowden and Clark J. Phys. C. vol. 14 p. L827 (1981)
c
c     By: Jens Jensen, Oersted Laboratory, Iniversity of Copenhagen
      REAL*8 A(3,3),B(3,3),TAU(16,3),R,PI,Q(3),QQ,U,V,UU,X,ER,
     /DD(2,16,16,3,3),QL(3),QS(3),W(3),VOL,QQL,QQS,SUM1,SUM2,UMAX0,
     /SUM(2,3,3)
      INTEGER I,N,J,NN,MM,K0,K1,K2,L0,L1,L2,LS,M0,M1,M2,LMAX,JJ,KK
      DATA PI/3.1415926536/
      READ(9,*)N
      WRITE(8,10)N
10    FORMAT(10X,'Number of sublattices:',I2,/)
      DO 11 I=1,N
11    READ(9,*)TAU(I,1),TAU(I,2),TAU(I,3)
c      TAU(1,1)=0.0
c      TAU(1,2)=0.0
c      TAU(1,3)=0.0
      DO 12 I=1,N
12    WRITE(8,13)I,TAU(I,1),TAU(I,2),TAU(I,3)
13    FORMAT(5X,'Tau(',I1,')',3(3X,F6.3))
      DO 14 I=1,3
14    READ(9,*)A(I,1),A(I,2),A(I,3)
      DO 15 I=1,3
15    WRITE(8,16)I,A(I,1),A(I,2),A(I,3)
16    FORMAT(7X,'A(',I1,')',3(3X,F6.3))
      DO 17 I=1,3
17    READ(9,*)B(I,1),B(I,2),B(I,3)
      DO 18 I=1,3
18    WRITE(8,19)I,B(I,1),B(I,2),B(I,3)
19    FORMAT(7X,'B(',I1,')',3(3X,F6.3))
      READ(9,*)R
      WRITE(8,20)R
20    FORMAT(/,10X,'Convergence factor:',F6.2,/)
      VOL=A(3,1)*(A(1,2)*A(2,3)-A(2,2)*A(1,3))+A(3,2)*(A(1,3)*A(2,1)-
     /A(1,1)*A(2,3))+A(3,3)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      WRITE(8,21)VOL
21    FORMAT(10X,'Volume of unit cell:',F7.3,/)
      READ(9,*)LMAX,UMAX0
      WRITE(8,23)LMAX,UMAX0
23    FORMAT(10X,'LMAX=',I3,3X,'UMAX0',F5.1,/)
100   READ(9,*)Q(1),Q(2),Q(3)
      U=0.0
      DO 22 I=1,3
22    U=U+Q(I)*Q(I)
      QQ=DSQRT(U)
      IF (QQ.LT.0.00001) GO TO 1000
      WRITE(8,24)Q(1),Q(2),Q(3)
24    FORMAT(/,/,10X,'Q = ',F7.3,3X,F7.3,3X,F7.3,/)
      DO 30 NN=1,N
      DO 30 MM=1,NN
      IF (NN.GT.1.AND.NN.EQ.MM) GO TO 30
      DO 40 J=1,3
      DO 40 K=1,J
      SUM1=0.0
      SUM2=0.0
      DO 50 K0=1,LMAX
      DO 50 K1=1,2
      K2=K0-1
      IF (K1.EQ.2) THEN
       IF (K2.EQ.0) GO TO 50
       K2=-K2
      END IF
      W(1)=FLOAT(K2)
      DO 51 L0=1,LMAX
      DO 51 L1=1,2
      L2=L0-1
      IF (L1.EQ.2) THEN
       IF (L2.EQ.0) GO TO 51
       L2=-L2
      END IF
      W(2)=FLOAT(L2)
      DO 52 M0=1,LMAX
      DO 52 M1=1,2
      M2=M0-1
      IF (M1.EQ.2) THEN
       IF (M2.EQ.0) GO TO 52
       M2=-M2
      END IF
      W(3)=FLOAT(M2)
      DO 60 JJ=1,3
      QL(JJ)=0.0
      DO 61 KK=1,3
61    QL(JJ)=QL(JJ)+W(KK)*B(KK,JJ)
      QS(JJ)=QL(JJ)+Q(JJ)
60    CONTINUE
      U=0.0
      DO 62 JJ=1,3
62    U=U+QS(JJ)*QS(JJ)
      QQS=DSQRT(U)
      U=QQS/(2.0*R)
      U=U*U
      IF (U.GT.UMAX0) GO TO 64
      UU=-4.0*PI*DEXP(-U)*QS(J)*QS(K)/(QQS*QQS)
      U=0.0
      IF (MM.NE.NN) THEN
       DO 63 JJ=1,3
63     U=U+QL(JJ)*(TAU(NN,JJ)-TAU(MM,JJ))
      END IF
      SUM1=SUM1+UU*DCOS(U)
      SUM2=SUM2-UU*DSIN(U)
64    CONTINUE
52    CONTINUE
51    CONTINUE
50    CONTINUE
      IF (J.EQ.K.AND.MM.EQ.NN) THEN
       DD(1,MM,MM,J,J)=4.0*R*R*R*VOL/(3*DSQRT(PI))+SUM1
       DD(2,MM,MM,J,J)=SUM2
       ELSE
       DD(1,MM,NN,J,K)=SUM1
       DD(2,MM,NN,J,K)=SUM2
       DD(1,MM,NN,K,J)=SUM1
       DD(2,MM,NN,K,J)=SUM2
      END IF
      SUM(1,J,K)=0.0
      SUM(2,J,K)=0.0
40    CONTINUE
      DO 70 K0=1,LMAX
      DO 70 K1=1,2
      K2=K0-1
      IF (K1.EQ.2) THEN
       IF (K2.EQ.0) GO TO 70
       K2=-K2
      END IF
      W(1)=FLOAT(K2)
      DO 71 L0=1,LMAX
      DO 71 L1=1,2
      L2=L0-1
      IF (L1.EQ.2) THEN
       IF (L2.EQ.0) GO TO 71
       L2=-L2
      END IF
      W(2)=FLOAT(L2)
      DO 72 M0=1,LMAX
      DO 72 M1=1,2
      M2=M0-1
      IF (M1.EQ.2) THEN
       IF (M2.EQ.0) GO TO 72
       M2=-M2
      END IF
      W(3)=FLOAT(M2)
      DO 80 JJ=1,3
      QL(JJ)=TAU(NN,JJ)-TAU(MM,JJ)
      DO 81 KK=1,3
81    QL(JJ)=QL(JJ)+W(KK)*A(KK,JJ)
80    CONTINUE
      U=0.0
      DO 82 JJ=1,3
82    U=U+QL(JJ)*QL(JJ)
      QQL=DSQRT(U)
      U=R*QQL
      U=U*U
      IF (U.LT.0.00001) GO TO 84
      IF (U.GT.UMAX0) GO TO 84
      UU=2.0*R*VOL*DEXP(-U)/(DSQRT(PI)*QQL*QQL)
      U=R*QQL
      V=1.0
      X=1.0
      LS=0
90    LS=LS+1
      V=V*(2.0*U*U)/(2.0*FLOAT(LS)+1.0)
      X=X+V
      IF (LS.LT.200) GO TO 90
      DO 41 J=1,3
      DO 41 K=1,J
      U=3.0*QL(J)*QL(K)/(QQL*QQL)
      IF (J.EQ.K) U=U-1.0
      ER=U*(UU*(1-X)+VOL/(QQL*QQL*QQL))+2.0*UU*R*R*QL(J)*QL(K)
      U=0.0
      DO 91 JJ=1,3
91    U=U+Q(JJ)*QL(JJ)
      SUM(1,J,K)=SUM(1,J,K)+ER*DCOS(U)
      SUM(2,J,K)=SUM(2,J,K)+ER*DSIN(U)
41    CONTINUE
84    CONTINUE
72    CONTINUE
71    CONTINUE
70    CONTINUE
      DO 42 J=1,3
      DO 42 K=1,J
       DD(1,MM,NN,J,K)=DD(1,MM,NN,J,K)+SUM(1,J,K)
       DD(2,MM,NN,J,K)=DD(2,MM,NN,J,K)+SUM(2,J,K)
       DD(1,MM,NN,K,J)=DD(1,MM,NN,J,K)
       DD(2,MM,NN,K,J)=DD(2,MM,NN,J,K)
42    CONTINUE
30    CONTINUE
      DO 92 NN=1,N
      DO 92 MM=1,NN
      IF (NN.GT.1.AND.NN.EQ.MM) GO TO 92
      WRITE(8,95)MM,NN
95    FORMAT(/,10X,'(',I2,';',I2,')')
       DO 93 J=1,3
93     WRITE(8,94)DD(1,MM,NN,J,1),DD(1,MM,NN,J,2),DD(1,MM,NN,J,3),
     /DD(2,MM,NN,J,1),DD(2,MM,NN,J,2),DD(2,MM,NN,J,3)
94    FORMAT(5X,3(2X,F10.4),10X,3(2X,F10.4))
92    CONTINUE
      GO TO 100
1000  STOP
      END
