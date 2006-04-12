C     Electron-electron bremsstrahlung, diff. cross section lab system
      DIMENSION S1(2),S2(2)
      REAL*8 A1,A2,ALPHA,ARC,AR0,CT,CT0,DTH,E1,E11,EN1,F1,F2,G0,G1,G2,
     1 G3,G4,G5,H1,H2,H3,H4,H5,H6,H7,H8,K,KMAX,KQ,L,L1,L2,L3,L4,MCQ,
     2 P1,P1Q,PEN,PI,PM,R1,R2,R44,RK,RO,RQ,RQ2,RQ4,RX1,RX2,RX3,S,S1,
     3 S2,TH,THO,THETA,W,W2,W4,W44,WIQ,WIQC,WQ,WQ2,WQ2Q,WQ4,WQ4Q,WR,
     4 WR1,WW,WWQ,X,X1,X1Q,X2,X2Q,X3,X3Q,X12,XQ,XR,ZPIA
      DATA ALPHA/0.7297353D-2/,AR0/0.5794675D-3/,MCQ/0.5109989D3/,
     1 PI/0.31415926536D1/,EN1/0.3D3/,PEN/0.1D3/,DTH/0.1D1/,NTH/180/
      NTH1=NTH+1
      ARC=PI/1.80D2
      ZPIA=0.2D1*PI*ALPHA
      WRITE(10,10) EN1
   10 FORMAT(' Electron-electron bremsstrahlung in the lab system for E1
     1 =',F6.1,' keV'/)
      E11=EN1/MCQ
      E1=E11+0.1D1
      P1Q=E11*(E11+0.2D1)
      P1=DSQRT(P1Q)
      WQ2=E1+E1
      WQ2Q=WQ2*WQ2
      WQ=WQ2+0.2D1
      W=DSQRT(WQ)
      WQ4=WQ2-0.2D1
      WQ4Q=WQ4*WQ4
      W44=0.25D0*WQ4
      WW=DSQRT(WQ4)
      WWQ=WQ*WQ4
      A1=ZPIA*WQ2/(W*WW)
      F1=(DEXP(A1)-0.1D1)/A1
      KMAX=E11/(E1-P1+0.1D1)
      PM=KMAX*MCQ
      WRITE(10,11) PM
   11 FORMAT(' Maximum photon energy =',F7.2,' keV'/)
      K=PEN/MCQ
      IF(K.GT.KMAX) GO TO 60
      KQ=K*K
      RK=0.1D1/K
      THO=180.D0
      CT0=(E1+0.1D1-E11/K)/P1
      IF(CT0.LT.-0.1D1) GO TO 15
      THO=DACOS(CT0)/ARC
   15 WRITE(10,16) PEN,THO
   16 FORMAT(' Photon energy =',F6.1,' keV, theta0 =',F6.1,' degrees:'/)
      THETA=0.D0
      DO 50 J=1,NTH1
      IF(THETA.GT.THO) GO TO 50
      TH=THETA*ARC
      CT=DCOS(TH)
      M=1
      X1=K*(E1-P1*CT)
      X2=K
      X1Q=X1*X1
      RX1=0.1D1/X1
      X3=X1
      X3Q=X1Q
      RX3=RX1
      X2Q=KQ
      RX2=RK
      X12=X1*X2
      X=X1+X2
      XQ=X*X
      H1=(X1-X2)/X
      RQ=WQ-X-X
      RO=DSQRT(RQ)
      RQ2=RQ-0.2D1
      RQ4=RQ-0.4D1
      R44=0.25D0*RQ4
      WR=DSQRT(RQ4)
      A2=ZPIA*RQ2/(RO*WR)
      F2=(DEXP(A2)-0.1D1)/A2
      XR=X12/RQ
      L1=DLOG(0.5D0*(RO+WR))
      W4=WW*DSQRT(W44*RQ4+0.4D1*XR)
      H2=W*WR+RO*WW
      L3=DLOG(0.125D0*H2*H2/X)
      L4=DLOG(0.1D1+(0.25D0*WR*W4+0.2D1*W44*R44)/XR)
      H3=0.125D0*RQ2/X*(RQ4+WQ)**2
      H4=0.4D1+0.8D1*WQ2/WQ4Q+0.5D0*RQ*RQ2/X+0.25D0*X/X12*(WQ2Q
     1+RQ2*RQ4-0.8D1*RQ2/WQ4)
      H5=0.2D1/WQ4Q*(WQ*WQ2*RQ4-RQ2-RQ2+0.4D1*RQ/WQ)
      H6=0.1D1+0.125D0*RQ2/X12*(WQ2Q+RQ2*RQ2-0.6D1*(WQ4+RQ)
     1+0.16D2*X/WQ4)+0.2D1/X12+RQ4/WQ4-0.8D1/WQ4Q
      H7=0.25D0*(RQ+WQ)/X12*H1*H1-0.25D0*(RX1-RX2)**2-0.5D0*RQ/XQ
     1+0.2D1/(WQ4Q*XR)*(0.1D1+WQ4)
      H8=RQ/WWQ*(0.12D2*WQ2Q/WQ4*X12-XQ*(0.8D1+0.48D2/WWQ))
   20 R1=RQ4+0.4D1*(X1+X1Q/RQ)
      WR1=DSQRT(R1)
      R2=RQ4+X1+X1
      W2=DSQRT(R44*X2Q+0.2D1*X*XR)
      L=DLOG(0.25D0*RO*RX1*(R2+WR*WR1))
      L2=DLOG(0.1D1+RQ*RX1/X*(R44*X2+0.5D0*WR*W2))
      G0=H7+H8/(X1Q*X1Q)-RQ*X/(W44*X1Q*X1)+RQ/X1Q*(0.4D1/WQ-0.15D1
     1-(0.8D1/WQ4-WQ2*X)/WWQ)-RQ*RX1/WQ4+((WQ-0.4D1*X2)/RQ
     2-0.25D0*WWQ/X12-0.4D1*RX1+(WQ-0.5D0*WQ2*RQ)/X1Q)/R1
      G1=0.4D1+(0.9D1-E1)*RX2-0.2D1*RX1*(WQ+WQ-X2+0.4D1)+0.75D0*WWQ/X12
     1+0.4D1*RX1/R2*(X2-0.3D1*X1+0.4D1+(RQ+RQ-0.6D1)*RX1)+RX1/R1
     2*(WQ4-0.4D1*XR)*(0.25D0*WQ*R2*RX2-RQ2-X1-X1)
      G2=RO*((RQ+0.2D1)/XQ+0.8D1/X1Q)
      G3=0.2D1*WQ2*RX2-(RQ2-X2)*RX1+0.5D0*RQ2*(RQ+X1)/X+WQ4/R2-0.2D1
     1+H3*RX1+RX1/R2*(RQ2*((WQ2-X)*X2-WQ2)-X2Q-X2Q-0.4D1*X2)+(RQ2*
     2(0.15D1*RQ4-0.5D0*WQ*(RQ-0.5D1))+X1Q+X1Q-0.6D1*X1+WQ2*X2)/(R2*X)
      G4=H4-0.15D1*WQ2*RX2+RX1*(X2Q-WQ-WQ-WQ*X2+X2+0.5D0*WQ2Q)
     1-0.5D0*X2*WQ2/X-0.25D0*RX2*(WQ4+RQ)/X*WQ2Q+0.4D1*XQ*RX1/WWQ
     2-(0.2D1*RQ2*X2+WQ4*X-0.4D1*WQ2+0.5D0*RX1*(0.8D1-RQ)*WQ2)/R2
     3+(0.5D0*RQ2*(0.3D1*RQ4-WQ*(RQ-0.5D1))-X1*(WQ-X1-X1+0.4D1))/(X*R2)
     4+H5*RX1+WQ2*X/(W44*X1Q)*(0.12D2*X/WWQ+X-E1)*(0.1D1-X*RX1/WQ)
      G5=H6-(0.2D1+X1+X1)*RX2
      S=WR*G0+G1*L/WR1+G2*L1+G3*L2/W2+0.2D1*G4*RO*L3/(W*WW*X1)
     1-G5*L4/W4
      S1(M)=S
      S2(M)=S/F2
      IF(M.EQ.2) GO TO 30
      M=2
      X1=X2
      X1Q=X2Q
      RX1=RX2
      X2=X3
      X2Q=X3Q
      RX2=RX3
      GO TO 20 
   30 WIQ=AR0*K*(S1(1)+S1(2))/(PI*RO*W*WW*MCQ)
      WIQC=AR0*K*F1*(S2(1)+S2(2))/(PI*RO*W*WW*MCQ)
      WRITE(10,40) THETA,WIQ,WIQC
   40 FORMAT(' theta =',F6.1,': CS =',E12.5,' b/(keV*sr), CSC =',
     1 E12.5,' b/(keV*sr)')
   50 THETA=THETA+DTH
      STOP
   60 WRITE(6,70)
   70 FORMAT(' The chosen photon energy is higher than the maximum photo
     1n energy')
      STOP  
      END 
