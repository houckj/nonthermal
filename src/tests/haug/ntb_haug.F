C
C  Date: Tue, 25 Nov 2003 10:08:31 +0100 (CET)
C  From: Eberhard Haug <haug@tat.physik.uni-tuebingen.de>
C  Subject: e-e bremsstrahlung
C  To: <houck@space.mit.edu>
C
C  Dr. Houck,
C
C  I do have software which evaluates the cross section of my 1975 paper,
C  both for the center-of-mass (cms) frame and the lab frame, and of the
C  integrated cross section of the A&A paper. In the attachment I send you
C  the FORTRAN codes (1) of the differential cross section in the cms and
C  (2) of the integrated cross section ion the cms.
C  The actual parameters of the codes reproduce (1) two curves of Fig. 2
C  and (2) a curve of Fig. 6 of my 1975 paper, Z. Naturforsch. 30a, 1099.
C  If you have further questions, please contact me.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c  J. Houck (Nov 2003)   code slightly modified and rearranged,
c                        moving most of the work into a subroutine
c
c electron-electron bremsstrahlung in the cms, photon spectrum
c       en = electron energy (keV)
c     phot = photon energy (keV)
c       wq = cross-section (b/keV)
c      wqc = cross-section w/ Coulomb correction (b/keV)
c
      double precision function eecms (en, phot)
      real*8 a1,a2,al,alp,al1,al2,al3,al4,al5,c,ct,dth,e,e1,e3,e4,
     & eek,ek,ekq,e2k,emk,en,ep,epp,epq,epq3,eq,fakt,f1,f2,fee,fk,g1,
     & g2,g3,g4,g5,g6,g7,g8,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,
     & h13,h14,h15,h16,k,kmax,kq,l,l2,l4,mcq,p,pck,phot,pi,pq,pqe,
     & pqr4,r1,r2,r4,r44,r44q,rek,ro,rp4,rq,s,si,ss,st,th,w2,w4,
     & we,wep,wp,wq,wqc,wr1,wr4,x1,x2,x12,y
      parameter (c=0.579468d-3,mcq=0.510999d3,fk=0.7297353d-2,
     &  pi=0.31415926536d1,nt=2000)
c     nt: number of theta values in numerical integration
      n1=nt-1
      dth=pi/nt
C      write(10,10) en
C 10   format('# electron-electron bremsstrahlung in the cms'/
C     &'# electron energy =',f7.1,' kev'/)
      e1=en/mcq
      e=e1+0.1d1
      pq=e1*(e1+0.2d1)
      p=dsqrt(pq)
      rp4=0.1d1/(pq*pq)
      eq=pq+0.1d1
      e3=e*eq
      e4=eq*eq
      ep=e*p
      epp=e+p
      epq=eq+pq
      y=epq+epq
      epq3=epq*epq*epq
      pqe=pq*e
      a1=pi*fk*epq/ep
      f1=(dexp(a1)-1)/a1
      al=dlog(epp)
      alp=al/ep
      kmax=pq/e
C
      k=phot/mcq
      if(k.ge.kmax) then
        eecms = 0.0
        return
      endif
C
      kq=k*k
      ek=e*k
      ekq=ek*ek
      emk=e-k
      rek=k/e
      h1=k/epp
      h2=k*epp
      h3=epq-ek-ek
      h4=pqe+k
      e2k=eq*k
      eek=eq-ek
      rq=0.4d1*eek
      ro=dsqrt(rq)
      we=0.5d0*ro
      r44=pq-ek
      wp=dsqrt(r44)
      wep=we*wp
      wr4=wp+wp
      r44q=r44*r44
      r4=0.4d1*r44
      pqr4=pq*r4
      a2=pi*fk*h3/wep
      f2=(dexp(a2)-0.1d1)/a2
      fee=f1/f2
      h5=(eq+0.1d1)/eek
      h6=ekq+0.5d0*ek-0.5d0-0.5d0*h3*(pqr4+eq*(kq+kq-0.2d1+ek/pq))
      h7=0.4d1*pq+0.16d2/(rek*rp4)-0.16d2+0.8d1*e4*h5
      h8=0.25d0*ek/pq+0.125d0*rp4-0.1d1
      h9=k*(e4+0.4d1*eq+0.25d1)/r44
      h10=pq*rq-kq
      h11=k*(0.3d1*e2k+0.8d1*e-k)/r44
      h12=0.2d1*epq3/e2k
      h13=r4+r4+ek+ek
      h14=ek+ek-0.4d1*eq-0.5d0+0.2d1*pqe/k
      h15=0.5d0/rek-0.5d0
      h16=pq-0.1d1-0.2d1*pq*ek-0.8d1*h15/rp4+(0.4d1*pqe-e+0.15d1/e)/k
      al1=dlog(we+wp)
      al2=dlog((p*ro+k)/(p*ro-k))
      al3=dlog((pqe+pqe+p*ro*wp)/k-epq)
      al4=dlog((ep+ep+h1)/(ep+ep-h2))
      al5=dlog((pq+pq-h1)/(pq+pq-h2))
      fakt=c/(pqe*phot)
      g1=0.4d1*e3+0.6d1*e+k+0.1d1/e+emk*((0.4d1-ek-ek)/pq+rp4+k/h4)
     &  +h9+h9
      g2=emk*(e/0.1875d0+(0.4d1*e-k-k)/pq+rp4/e)+0.4d1*kq+rek
     &  +(eq+eq+0.2d1)/r44
      g3=0.5d0*ek/h4-0.4d1*e+e2k*h5-h9
      g4=pq*e2k/h4*(e3+e3+k)-0.3d1*e2k-0.5d1*k-0.4d1*pqe*(eq+0.1d1)
      g5=0.18d2-rek-rek+0.1d1/eq+(ek*(0.8d1*e4-0.4d1*eq-0.28d2)
     &  +0.12d2*kq+0.16d2*eq*(0.2d1-e4)+h11+h11)/h10
      g6=(0.32d2*pq*emk-0.19d2*k+0.14d2/e)/0.3d1+0.8d1*e*(kq-0.1d1)
     &  +k*rek-k/eq-k/pqe*(emk+0.5d0/pqe)-epq/(e*r44)+(0.4d1*e-e2k-e2k
     &  +k/h10*(0.8d1*eq-0.12d2*ek*eek-(e4+e4+pq)*kq+h11))/eek
     &  +alp*(0.16d2*pqe-0.10d2*pq*k-0.4d1*e-k+0.55d1/e+(0.4d1*k-e)/pq
     &  +k*rp4+k*(0.15d1*epq-0.2d1*eq*ek)/r44+0.5d0*pq*(e+k)/r44q-h12)
      g6=g6+al1/wep*((0.18d2*k-0.28d2*e)*pq-0.6d1*e*kq-0.8d1*emk
     &  -0.4d1/e+h12-(e2k+k)/r44)+al5/p*(eq+eq-0.4d1+rek+0.75d0/eq
     &  +(ekq-0.75d0*rek*epq-0.25d0*pq*(0.1d1+rek)/r44)/r44)
      s=(g1*al-g2*p)*wp/we+g3*wp*al2+0.5d0*g4*al4/wep+g5*p*al1+g6*al3
      th=dth
      si=0.4d1
      ss=0.d0
      do n=1,n1
        ct=dcos(th)
        st=dsin(th)
        pck=p*ct*k
        x1=ek-pck
        x2=ek+pck
        x12=x1*x2
        r1=0.4d1*(pq-pck)+x1*x1/eek
        wr1=dsqrt(r1)
        r2=pq+pq-x2
        w2=dsqrt(x2*(r44*x2+k*x1/emk))
        w4=p*dsqrt(pqr4+x12/eek)
        l=dlog(we/x1*(r2+wp*wr1))
        l2=dlog(0.1d1+emk/(k*x1)*(r44*(x2+x2)+wr4*w2))
        l4=dlog(0.1d1+(eek+eek)/x12*(pqr4+wr4*w4))
        g7=l/wr1*((h13+0.2d1/x1)/(x1*r2)+h7/(x2*r1))
        g8=l2/w2*(h14+h15*x1+y/x2+(h16-epq/x1)/r2)+l4/w4*(h8+h6/x12)
        ss=ss+(g7+g8)*si*st
        if(si.eq.0.4d1) then
          si=0.2d1
        else
          si=0.4d1
        endif
        th=th+dth
      enddo
      ss=ss*p*kq*dth/(0.3d1*we)
      wq=fakt*(s+ss)
      wqc=fee*wq

C      write(10,70) phot,wq,wqc
C  70  format (1x, e12.5, 1x, e12.5, 1x, e12.5)
c 70   format(' photon energy =',f6.1,' kev: cs =',e12.5,' b/kev',
c     & 3x,'csc =',e12.5,' b/kev')
c     cs: cross section
c     csc: cross section including coulomb correction

      eecms = wqc

      return
      end

#ifdef TESTING
      program ntbrems
      double precision elg, xlg, phot, en, cms, eecms

      do elg = 0.0d0, 6.0d0, 0.025d0
        en = 10.0**elg
        do xlg=0.0d0, elg+0.5, 0.025d0
          phot = 10.0d0**xlg
          cms = eecms (en, phot)
          write(10,70) en, phot, cms
  70      format (1x, 1pe12.5, 1x, 1pe12.5, 1x, 1pe12.5)
        enddo
      enddo

      stop
      end
#endif
