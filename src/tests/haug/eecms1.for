C     Electron-electron bremsstrahlung in the cms, photon spectrum
      REAL*8 A1,A2,AL,ALP,AL1,AL2,AL3,AL4,AL5,C,CT,DP,DTH,E,E1,E3,E4,
     1 EEK,EK,EKQ,E2K,EMK,EN,EP,EPP,EPQ,EPQ3,EQ,FAKT,F1,F2,FEE,FK,G1,
     2 G2,G3,G4,G5,G6,G7,G8,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12,
     3 H13,H14,H15,H16,K,KMAX,KQ,L,L2,L4,MCQ,P,PCK,PHOT,PI,PQ,PQE,
     4 PQR4,R1,R2,R4,R44,R44Q,REK,RO,RP4,RQ,S,SI,SS,ST,TH,W2,W4,
     5 WE,WEP,WP,WQ,WQC,WR1,WR4,X1,X2,X12,Y
      DATA C/0.579468D-3/,MCQ/0.510999D3/,FK/0.7297353D-2/,NP/200/,
     1 PI/0.31415926536D1/,EN/0.100D3/,NT/2000/
C     EN: electron energy (keV)
C     NT: number of theta values in numerical integration 
C     NP: number of photon energies
      N1=NT-1
      DP=0.1D-1*EN
      DTH=PI/NT
      WRITE(10,10) EN
 10   FORMAT(' Electron-electron bremsstrahlung in the cms'/
     1' electron energy =',F7.1,' keV'/)
      E1=EN/MCQ
      E=E1+0.1D1
      PQ=E1*(E1+0.2D1)
      P=DSQRT(PQ)
      RP4=0.1D1/(PQ*PQ)
      EQ=PQ+0.1D1
      E3=E*EQ
      E4=EQ*EQ
      EP=E*P
      EPP=E+P
      EPQ=EQ+PQ
      Y=EPQ+EPQ
      EPQ3=EPQ*EPQ*EPQ
      PQE=PQ*E
      A1=PI*FK*EPQ/EP
      F1=(DEXP(A1)-1)/A1
      AL=DLOG(EPP)
      ALP=AL/EP
      KMAX=PQ/E
      PHOT=DP
      DO 80 JP=1,NP
      K=PHOT/MCQ
      IF(K.GE.KMAX) GO TO 80
      KQ=K*K
      EK=E*K
      EKQ=EK*EK
      EMK=E-K
      REK=K/E
      H1=K/EPP
      H2=K*EPP
      H3=EPQ-EK-EK
      H4=PQE+K
      E2K=EQ*K
      EEK=EQ-EK
      RQ=0.4D1*EEK
      RO=DSQRT(RQ)
      WE=0.5D0*RO
      R44=PQ-EK
      WP=DSQRT(R44)
      WEP=WE*WP
      WR4=WP+WP
      R44Q=R44*R44
      R4=0.4D1*R44
      PQR4=PQ*R4
      A2=PI*FK*H3/WEP
      F2=(DEXP(A2)-0.1D1)/A2
      FEE=F1/F2
      H5=(EQ+0.1D1)/EEK
      H6=EKQ+0.5D0*EK-0.5D0-0.5D0*H3*(PQR4+EQ*(KQ+KQ-0.2D1+EK/PQ))
      H7=0.4D1*PQ+0.16D2/(REK*RP4)-0.16D2+0.8D1*E4*H5
      H8=0.25D0*EK/PQ+0.125D0*RP4-0.1D1
      H9=K*(E4+0.4D1*EQ+0.25D1)/R44
      H10=PQ*RQ-KQ
      H11=K*(0.3D1*E2K+0.8D1*E-K)/R44
      H12=0.2D1*EPQ3/E2K
      H13=R4+R4+EK+EK
      H14=EK+EK-0.4D1*EQ-0.5D0+0.2D1*PQE/K
      H15=0.5D0/REK-0.5D0
      H16=PQ-0.1D1-0.2D1*PQ*EK-0.8D1*H15/RP4+(0.4D1*PQE-E+0.15D1/E)/K
      AL1=DLOG(WE+WP)
      AL2=DLOG((P*RO+K)/(P*RO-K))
      AL3=DLOG((PQE+PQE+P*RO*WP)/K-EPQ)
      AL4=DLOG((EP+EP+H1)/(EP+EP-H2))
      AL5=DLOG((PQ+PQ-H1)/(PQ+PQ-H2))
      FAKT=C/(PQE*PHOT)
      G1=0.4D1*E3+0.6D1*E+K+0.1D1/E+EMK*((0.4D1-EK-EK)/PQ+RP4+K/H4)
     1+H9+H9
      G2=EMK*(E/0.1875D0+(0.4D1*E-K-K)/PQ+RP4/E)+0.4D1*KQ+REK
     1+(EQ+EQ+0.2D1)/R44
      G3=0.5D0*EK/H4-0.4D1*E+E2K*H5-H9
      G4=PQ*E2K/H4*(E3+E3+K)-0.3D1*E2K-0.5D1*K-0.4D1*PQE*(EQ+0.1D1)
      G5=0.18D2-REK-REK+0.1D1/EQ+(EK*(0.8D1*E4-0.4D1*EQ-0.28D2)
     1+0.12D2*KQ+0.16D2*EQ*(0.2D1-E4)+H11+H11)/H10
      G6=(0.32D2*PQ*EMK-0.19D2*K+0.14D2/E)/0.3D1+0.8D1*E*(KQ-0.1D1)
     1+K*REK-K/EQ-K/PQE*(EMK+0.5D0/PQE)-EPQ/(E*R44)+(0.4D1*E-E2K-E2K
     2+K/H10*(0.8D1*EQ-0.12D2*EK*EEK-(E4+E4+PQ)*KQ+H11))/EEK
     3+ALP*(0.16D2*PQE-0.10D2*PQ*K-0.4D1*E-K+0.55D1/E+(0.4D1*K-E)/PQ
     4+K*RP4+K*(0.15D1*EPQ-0.2D1*EQ*EK)/R44+0.5D0*PQ*(E+K)/R44Q-H12)
      G6=G6+AL1/WEP*((0.18D2*K-0.28D2*E)*PQ-0.6D1*E*KQ-0.8D1*EMK
     1-0.4D1/E+H12-(E2K+K)/R44)+AL5/P*(EQ+EQ-0.4D1+REK+0.75D0/EQ
     2+(EKQ-0.75D0*REK*EPQ-0.25D0*PQ*(0.1D1+REK)/R44)/R44)
      S=(G1*AL-G2*P)*WP/WE+G3*WP*AL2+0.5D0*G4*AL4/WEP+G5*P*AL1+G6*AL3
      TH=DTH
      SI=0.4D1
      SS=0.
      DO 50 N=1,N1
      CT=DCOS(TH)
      ST=DSIN(TH)
      PCK=P*CT*K
      X1=EK-PCK
      X2=EK+PCK
      X12=X1*X2
      R1=0.4D1*(PQ-PCK)+X1*X1/EEK
      WR1=DSQRT(R1)
      R2=PQ+PQ-X2
      W2=DSQRT(X2*(R44*X2+K*X1/EMK))
      W4=P*DSQRT(PQR4+X12/EEK)
      L=DLOG(WE/X1*(R2+WP*WR1))
      L2=DLOG(0.1D1+EMK/(K*X1)*(R44*(X2+X2)+WR4*W2))
      L4=DLOG(0.1D1+(EEK+EEK)/X12*(PQR4+WR4*W4))
      G7=L/WR1*((H13+0.2D1/X1)/(X1*R2)+H7/(X2*R1))
      G8=L2/W2*(H14+H15*X1+Y/X2+(H16-EPQ/X1)/R2)+L4/W4*(H8+H6/X12)
      SS=SS+(G7+G8)*SI*ST
      IF(SI.EQ.0.4D1) GO TO 40
      SI=0.4D1
      GO TO 50
 40   SI=0.2D1
 50   TH=TH+DTH
      SS=SS*P*KQ*DTH/(0.3D1*WE)
      WQ=FAKT*(S+SS)
      WQC=FEE*WQ
      WRITE(10,70) PHOT,WQ,WQC
 70   FORMAT(' Photon energy =',F6.1,' keV: CS =',E12.5,' b/keV',
     1 3X,'CSC =',E12.5,' b/keV')
C     CS: cross section
C     CSC: cross section including Coulomb correction
      IF(JP.EQ.10) DP=0.5D1*DP
      IF(JP.EQ.18) DP=0.2D1*DP
 80   PHOT=PHOT+DP
      STOP
      END
