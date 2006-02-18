C-*-fortran-*-
c
c John C. Houck <houck@space.mit.edu> made a few changes,
c all of which are marked with JCH

      integer function ibsearch (t, x, n)
      double precision t, xt, x(n)
      integer n0,n1,n2

      if ((t.lt.x(1)) .or. (x(n).le.t)) then
        ibsearch = 0
        return
      endif

      n0 = 1
      n1 = n

      do while (n1.gt.n0+1)
        n2 = (n0 + n1) / 2
        xt = x(n2)
        if (t .le. xt) then
          if (xt .eq. t) then
            ibsearch = n2
            return
          endif
          n1 = n2
        else
          n0 = n2
        endif
      enddo

      ibsearch = n0
      
      return
      end
c
c   VERSION 2.2
c
c   This library contains routines for B-spline interpolation in
c   one, two, and three dimensions. Part of the routines are based
c   on the book by Carl de Boor: A practical guide to Splines (Springer,
c   New-York 1978) and have the same calling sequence and names as
c   the corresponding routines from the IMSL library. For documen-
c   tation see the additional files. NOTE: The results in the demo
c   routines may vary slightly on different architectures.
c
c   by W. Schadow 07/19/98
c   last changed by W. Schadow 07/28/2000
c
c
c   Wolfgang Schadow
c   TRIUMF
c   4004 Wesbrook Mall
c   Vancouver, B.C. V6T 2A3
c   Canada
c
c   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
c
c   www  : http://www.triumf.ca/people/schadow
c
c
c  ------------------------------------------------------------------
c
c
c   Copyright (C) 2000 Wolfgang Schadow
c
c   This library is free software; you can redistribute it and/or
c   modify it under the terms of the GNU Library General Public
c   License as published by the Free Software Foundation; either
c   version 2 of the License, or (at your option) any later version.
c
c   This library is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
c   Library General Public License for more details.
c
c   You should have received a copy of the GNU Library General Public
c   License along with this library; if not, write to the
c   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
c   Boston, MA  02111-1307, USA.
c
c
c  ------------------------------------------------------------------
c
c
c   The following routines are included:
c
c            dbsnak
c
c            dbsint
c            dbsval
c            dbsder
c            dbs1gd
c
c            dbs2in
c            dbs2dr
c            dbs2vl
c            dbs2gd
c
c            dbs3in
c            dbs3vl
c            dbs3dr
c            dbs3gd
c
c  ------------------------------------------------------------------
c
c  NEW: corrected some error messages
c       some changes in the checks of dbs3dg to find a possible
c       error earlier.
c
c  ------------------------------------------------------------------
c
c  NEW: documentation included, changed some comments
c
c  ------------------------------------------------------------------
c

c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbsnak(nx,xvec,kxord,xknot)

c
c  Compute the `not-a-knot' spline knot sequence.
c  (see de Boor p. 167)
c
c   nx     - number of data points.  (input)
c   xvec   - array of length ndata containing the location of the
c            data points.  (input)
c   kxord  - order of the spline.  (input)
c   xknot  - array of length ndata+korder containing the knot
c            sequence.  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xvec(nx),xknot(nx+kxord)

      logical first

      save first,eps

      data first/.true./

      if (first) then
         first=.false.
         eps = dlamch('precision')
cJCH         write(6,*) 'subroutine dbsnak: '
cJCH         write(6,*) 'eps = ',eps
      endif

      if((kxord .lt. 0) .or. (kxord .gt. nx)) then
         write(6,*) 'subroutine dbsnak: error'
         write(6,*) '0 <= kxord <= nx is required.'
         write(6,*) 'kxord = ', kxord, ' and nx = ', nx,  ' is given.'
         stop
      endif

      do 30 i = 1, kxord
         xknot(i) = xvec(1)
 30   continue

      if(mod(kxord,2) .eq. 0) then
         do 40 i = kxord+1,nx
            xknot(i) = xvec(i-kxord/2)
 40      continue
      else
         do 50 i = kxord+1,nx
            xknot(i) = 0.5d0 * (xvec(i-kxord/2) + xvec(i-kxord/2-1))
 50      continue
      endif

c      do 60 i = nx+1,nx+kxord
c         xknot(i) = xvec(nx) * (1.0d0 + eps)
c 60   continue
cJCH  original algorithm fails when grid points are not positive
      do 60 i = nx+1,nx+kxord
        if (xvec(nx) .ne. 0.0d0) then
          xknot(i) = xvec(nx) + eps * abs(xvec(nx))
        else
          xknot(i) = xvec(nx) + eps * abs(xvec(nx-1))
        endif
 60   continue

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef)

c
c  Computes the spline interpolant, returning the B-spline coefficients.
c  (see de Boor p. 204)
c
c   nx     - number of data points.  (input)
c   xvec   - array of length nx containing the data point
c            abscissas.  (input)
c   xdata  - array of length ndata containing the data point
c            ordinates.  (input)
c   kx     - order of the spline.  (input)
c            korder must be less than or equal to ndata.
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   bscoef - array of length ndata containing the B-spline
c            coefficients.  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,nxmax=200)

      dimension xdata(nx),xvec(nx),xknot(nx+kx),bcoef(nx)
      dimension work((2*kxmax-1)*nxmax)

      if ((kx .gt. kxmax) .or. (nx .gt. nxmax)) then
         write(6,*) 'subroutine dbsint: error'
         write(6,*) 'kx > kxmax or nx > nxmax'
         write(6,*) 'kx = ',kx,'  kxmax = ',kxmax
         write(6,*) 'nx = ',nx,'  nxmax = ',nxmax
         stop
      endif

      nxp1  = nx + 1
      kxm1  = kx - 1
      kpkm2 = 2 * kxm1
      leftx = kx
      lenq  = nx * (kx + kxm1)

      do 10 i = 1, lenq
         work(i) = 0.d0
 10   continue

      do 20 ix = 1,nx
         xveci  = xvec(ix)
         ilp1mx = min0(ix+kx,nxp1)
         leftx   = max0(leftx,ix)
         if (xveci .lt. xknot(leftx)) goto 998
 30      if (xveci .lt. xknot(leftx+1)) go to 40
         leftx = leftx + 1
         if (leftx .lt. ilp1mx) go to 30
         leftx = leftx - 1
         if (xveci .gt. xknot(leftx+1)) goto 998
 40      call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
         jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
         do 50 ik = 1,kx
            jj       = jj + kpkm2
            work(jj) = bcoef(ik)
 50      continue
 20   continue

      call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)
      go to (60,999), iflag

 60   do 70 ix = 1,nx
         bcoef(ix) = xdata(ix)
 70   continue

      call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)

      return

 998  write(6,*) 'subroutine dbsint:'
      write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
      write(6,*) ix,xknot(ix),xknot(ix+1)

      stop

 999  write(6,*) 'subroutine dbsint: error'
      write(6,*) 'no solution of linear equation system !!!'

      stop

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsval(x,kx,xknot,nx,bcoef)

c
c  Evaluates a spline, given its B-spline representation.
c
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   dbsval - value of the spline at x.  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax)


c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbsval:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

c
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c

C       leftx = 0

c      do 10 ix = 1,nx+kx-1
c         if (xknot(ix) .gt. xknot(ix+1)) then
c             write(6,*) 'subroutine dbsval:'
c             write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
c             write(6,*) ix,xknot(ix),xknot(ix+1)
c          stop
c          endif
c         if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
c10    continue

cJCH  replaced linear search with binary search
      leftx = ibsearch (x, xknot, nx+kx)

      if(leftx .eq. 0) then
         write(6,*) 'subroutine dbsval:'
         write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
         write(6,*) 'x = ', x
         stop
      endif

      do 20 ik = 1,kx-1
         work(ik) = bcoef(leftx+ik-kx)
         dl(ik)   = x - xknot(leftx+ik-kx)
         dr(ik)   = xknot(leftx+ik) - x
 20   continue

      work(kx)  = bcoef(leftx)
      dl(kx)    = x - xknot(leftx)

      do 30 ik = 1,kx-1
         save2 = work(ik)
         do 40  il = ik+1,kx
            save1 = work(il)
            work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .           / (dl(il) + dr(il - ik))
            save2 = save1
 40      continue
 30   continue

      dbsval = work(kx)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsder(iderx,x,kx,xknot,nx,bcoef)

c
c  Evaluates the derivative of a spline, given its B-spline representation.
c
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   dbsder - value of the iderx-th derivative of the spline at x.
c            (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax),bsp(kxmax)

c
c     check if <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbsder:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

c
c     check if xknot(i) <= xknot(i+1) and calculation of i so that
c     xknot(i) <= x < xknot(i+1)
c

      leftx = 0
      do 10 ix = 1,nx+kx-1
         if (xknot(ix) .gt. xknot(ix+1)) then
            write(6,*) 'subroutine dbsder:'
            write(6,*) 'xknot(ix) <= xknot(ix+1) required.'
            stop
         endif
         if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
 10   continue

      if (leftx .eq. 0) then
         write(6,*) 'subroutine dbsder:'
         write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
         write(6,*) 'xknot(1)     = ', xknot(1)
         write(6,*) 'xknot(nx+kx) = ', xknot(nx+kx)
         write(6,*) '         x   = ', x
         stop
      endif

      if (iderx .eq. 0) then

         do 20 ik = 1,kx-1
            work(ik) = bcoef(leftx+ik-kx)
            dl(ik)   = x - xknot(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
 20      continue

         work(kx)  = bcoef(leftx)
         dl(kx)    = x - xknot(leftx)

         do 30 ik = 1,kx-1
            save2 = work(ik)
            do 40  il = ik+1,kx
               save1 = work(il)
               work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .              / (dl(il) + dr(il - ik))
               save2 = save1
 40         continue
 30      continue

         dbsder = work(kx)

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
         bsp(1) = 1.0d0
         do 50 ik = 1,kx-iderx-1
            dr(ik) = xknot(leftx+ik) - x
            dl(ik) = x - xknot(leftx+1-ik)
            save   = bsp(1)
            bsp(1) = 0.0d0
            do 60 il = 1,ik
               y         = save / (dr(il) + dl(ik+1-il))
               bsp(il)   = bsp(il) + dr(il) * y
               save      = bsp(il+1)
               bsp(il+1) = dl(ik+1-il) * y
 60         continue
 50      continue

         do 70 ik = 1,kx
            work(ik) = bcoef(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
            dl(ik)   = x - xknot(leftx+ik-kx)
 70      continue

         do 80 ik = 1,iderx
            dik   = dble(kx - ik)
            save2 = work(ik)
            do 90  il = ik+1,kx
               save1    = work(il)
               work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
               save2    = save1
 90         continue
 80      continue

         sum = 0.0d0

         do 100 i = 1,kx-iderx
            sum = sum + bsp(i) * work(iderx+i)
 100     continue

         dbsder = sum

      else
         dbsder = 0.0d0
      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs1gd(iderx,nxvec,xvec,kx,xknot,nx,bcoef,val)

c
c  Evaluates the derivative of a spline on a grid, given its B-spline
c  representation.
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   nxvec  - length of vector xvec.  (input)
c   xvec   - array of length nxvec containing the points at which the
c            spline is to be evaluated.  (input)
c            xvec should be strictly increasing.
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   val    - array of length nxvec containing the values of the
c            iderx-th derivative of the spline at the points in
c            xvec.  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,nxmax=200)

      dimension xvec(nxvec),xknot(nx+kx),bcoef(nx),val(nxvec)
      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),term(nxmax),save2(nxmax)
      dimension leftx(nxmax),work(nxmax,kxmax)

      logical same,next

c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs1gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 10 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 10   continue

      do 20 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs1gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 20   continue

      do 30 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs1gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 30   continue

      if (iderx .eq. 0) then

         do 40 ix = 1,nxvec
            biatx(ix,1) = 1.d0
            val(ix)     = 0.d0
 40      continue

         do 50 ik = 1,kx-1
            do 60 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 60         continue

            do 70 il = 1,ik
               do 80 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 80            continue
 70         continue

            do 90 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 90         continue
 50      continue

         do 100 ik = 1,kx
            do 110 ix = 1,nxvec
               val(ix) = val(ix) + biatx(ix,ik) * bcoef(leftx(ix)-kx+ik)
 110        continue
 100     continue

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

         do 120 ix = 1,nxvec
            biatx(ix,1) = 1.d0
            val(ix)     = 0.d0
 120     continue

         do 130 ik = 1,kx-iderx-1
            do 140 ix = 1,nxvec
               dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix)    = biatx(ix,1)
               biatx(ix,1) = 0.0d0
               do 150 il = 1,ik
                  term(ix)       = save1(ix)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il)   = biatx(ix,il) + dr(ix,il) * term(ix)
                  save1(ix)      = biatx(ix,il+1)
                  biatx(ix,il+1) = dl(ix,ik+1-il) * term(ix)
 150           continue
 140        continue
 130     continue

         do 160 ik = 1,kx
            do 170 ix = 1,nxvec
               work(ix,ik) = bcoef(leftx(ix)+ik-kx)
               dr(ix,ik)   = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik)   = xvec(ix) - xknot(leftx(ix)+ik-kx)
 170        continue
 160     continue

         do 180 ik = 1,iderx
            dik   = dble(kx - ik)
            do 190 ix = 1,nxvec
               save2(ix) = work(ix,ik)
               do 200  il = ik+1,kx
                  save1(ix)   = work(ix,il)
                  work(ix,il) = dik * (work(ix,il) - save2(ix))
     .                 /(dl(ix,il) + dr(ix,il-ik))
                  save2(ix)   = save1(ix)
 200           continue
 190        continue
 180     continue

         do 210 i = 1,kx-iderx
            do 220 ix = 1,nxvec
               val(ix) = val(ix) + biatx(ix,i) * work(ix,iderx+i)
 220        continue
 210     continue

      else

         do 230 ix = 1,nxvec
            val(ix) = 0.0d0
 230     continue

      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbsdca(iderx,x,kx,xknot,nx,bcoef,leftx)

c
c This routine is equivalent to the routine dbsder, but it does not
c check the parameters!!c
c
c Evaluates the derivative of a spline, given its B-spline representation.
c
c
c   iderx  - order of the derivative to be evaluated.  (input)
c            in particular, iderx = 0 returns the value of the
c            spline.
c   x      - point at which the spline is to be evaluated.  (input)
c   kx     - order of the spline.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence.  (input)
c            xknot must be nondecreasing.
c   nx     - number of B-spline coefficients.  (input)
c   bcoef  - array of length nx containing the B-spline
c            coefficients.  (input)
c   leftx  - number of the intervall of xknot that includes x
c   dbsdca - value of the ideriv-th derivative of the spline at x.
c            (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8)

      dimension xknot(nx+kx),bcoef(nx)
      dimension work(kxmax),dl(kxmax),dr(kxmax),bsp(kxmax)


      if (iderx .eq. 0) then

         do 20 ik = 1,kx-1
            work(ik) = bcoef(leftx+ik-kx)
            dl(ik)   = x - xknot(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
 20      continue

         work(kx)  = bcoef(leftx)
         dl(kx)    = x - xknot(leftx)

         do 30 ik = 1,kx-1
            save2 = work(ik)
            do 40  il = ik+1,kx
               save1 = work(il)
               work(il) = (dl(il) * work(il) + dr(il-ik) * save2)
     .              / (dl(il) + dr(il - ik))
               save2 = save1
 40         continue
 30      continue

         dbsdca = work(kx)

      elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then
         bsp(1) = 1.0d0
         do 50 ik = 1,kx-iderx-1
            dr(ik) = xknot(leftx+ik) - x
            dl(ik) = x - xknot(leftx+1-ik)
            save   = bsp(1)
            bsp(1) = 0.0d0
            do 60 il = 1,ik
               y         = save / (dr(il) + dl(ik+1-il))
               bsp(il)   = bsp(il) + dr(il) * y
               save      = bsp(il+1)
               bsp(il+1) = dl(ik+1-il) * y
 60         continue
 50      continue

         do 70 ik = 1,kx
            work(ik) = bcoef(leftx+ik-kx)
            dr(ik)   = xknot(leftx+ik) - x
            dl(ik)   = x - xknot(leftx+ik-kx)
 70      continue

         do 80 ik = 1,iderx
            dik   = dble(kx - ik)
            save2 = work(ik)
            do 90  il = ik+1,kx
               save1    = work(il)
               work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
               save2    = save1
 90         continue
 80      continue

         sum = 0.0d0

         do 100 i = 1,kx-iderx
            sum = sum + bsp(i) * work(iderx+i)
 100     continue

         dbsdca = sum

      else
         dbsdca = 0.0d0
      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs2in(nx,xvec,ny,yvec,xydata,ldf,kx,ky,xknot,
     .     yknot,bcoef)

c
c  Computes a two-dimensional tensor-product spline interpolant,
c  returning the tensor-product B-spline coefficients.
c
c    nx     - number of data points in the x-direction.  (input)
c    xvec   - array of length nx containing the data points in
c             the x-direction.  (input)
c             xdata must be strictly increasing.
c    ny     - number of data points in the y-direction.  (input)
c    yvec   - array of length ny containing the data points in
c             the y-direction.  (input)
c             ydata must be strictly increasing.
c    xydata - array of size nx by nydata containing the values to
c             be interpolated.  (input)
c             fdata(i,j) is the value at (xdata(i),ydata(j)).
c    ldf    - the leading dimension of fdata exactly as specified in
c             the dimension statement of the calling program.
c             (input)
c    kx     - order of the spline in the x-direction.  (input)
c             kxord must be less than or equal to nxdata.
c    ky     - order of the spline in the y-direction.  (input)
c             kyord must be less than or equal to nydata.
c    xknot  - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c    yknot  - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c    bcoef  - array of length nx*ny containing the
c             tensor-product B-spline coefficients.  (output)
c             bscoef is treated internally as a matrix of size nxdata
c             by nydata.
c

      implicit double precision (a-h,o-z), integer (i-n)

c$$$      parameter(nxmax=130,kxmax=8,nymax=130,kymax=8)
c FIXME      
      parameter(nxmax=1024,kxmax=8,nymax=1024,kymax=8)

c
c     dimensions should be
c                  work1(max(nx,ny),max(nx,ny))
c                  work2(max(nx,ny))
c                  work3(max((2*kx-1)*nx,(2*ky-1)*ny))
c

      dimension work1(nxmax,nymax),work2(nxmax)
      dimension work3((2*kxmax-1)*nxmax)

      dimension xvec(nx),xknot(nx+kx),yvec(ny),yknot(ny+ky)
      dimension xydata(ldf,*),bcoef(nx,ny)


      if ((kx .gt. kxmax) .or. (nx .gt. nxmax)) then
         write(6,*) 'subroutine dbs2in: error'
         write(6,*) 'kx > kxmax or nx > nxmax'
         write(6,*) 'kx = ',kx,'  kxmax = ',kxmax
         write(6,*) 'nx = ',nx,'  nxmax = ',nxmax
         stop
      endif

      if ((ky .gt. kymax) .or. (ny .gt. nymax)) then
         write(6,*) 'subroutine dbs2in: error'
         write(6,*) 'ky > kymax or ny > nymax'
         write(6,*) 'ky = ',ky,'  kymax = ',kymax
         write(6,*) 'ny = ',ny,'  nymax = ',nymax
         stop
      endif

      call spli2d(xvec,ldf,xydata,xknot,nx,kx,ny,work2,work3,work1)
      call spli2d(yvec,ny, work1, yknot,ny,ky,nx,work2,work3,bcoef)

      return
      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine spli2d(xyvec,ld,xydata,xyknot,n,k,m,work2,work3,bcoef)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xyvec(n),xyknot(n+k),xydata(ld,m),bcoef(m,n)
      dimension work2(n),work3((2*k-1)*n)

      np1   = n + 1
      km1   = k - 1
      kpkm2 = 2 * km1
      left  = k
      lenq  = n * (k + km1)

      do 10 i = 1,lenq
         work3(i) = 0.0d0
 10   continue

      do 20 i = 1,n
         xyveci  = xyvec(i)
         ilp1mx = min0(i+k,np1)
         left   = max0(left,i)
         if (xyveci .lt. xyknot(left)) go to 998
 30      if (xyveci .lt. xyknot(left+1)) go to 40
         left = left + 1
         if (left .lt. ilp1mx) go to 30
         left = left - 1
         if (xyveci .gt. xyknot(left+1)) go to 998
 40      call bsplvb(xyknot,n+k,k,1,xyveci,left,work2)
         jj = i - left + 1 + (left - k) * (k + km1)
         do 50 j = 1,k
            jj        = jj + kpkm2
            work3(jj) = work2(j)
 50      continue
 20   continue

      call banfac(work3,k+km1,n,km1,km1,iflag )

      go to (60,999), iflag

 60   do 70 j = 1,m
         do 80 i = 1,n
            work2(i) = xydata(i,j)
 80      continue

         call banslv(work3,k+km1,n,km1,km1,work2)

         do 90 i = 1,n
            bcoef(j,i) = work2(i)
 90      continue
 70   continue

      return

 998  write(6,*) 'subroutine db2in:'
      write(6,*) 'i with knot(i) <= x/y < knot(i+1) required.'
      write(6,*) 'knot(1)   = ', xyknot(1)
      write(6,*) 'knot(n+k) = ', xyknot(n+k)
      write(6,*) '      x/y = ', xyveci

      stop

 999  write(6,*) 'subroutine dbs2in: error'
      write(6,*) 'no solution of linear equation system !!!'

      stop

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

c
c  evaluates a two-dimensional tensor-product spline, given its
c  tensor-product B-spline representation.    use numeric
c
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   bcoef  - array of length nx*ny containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny.
c   dbs2vl - value of the spline at (x,y).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)

      dimension xknot(nx+kx),yknot(ny+ky),bcoef(nx,ny)
      dimension work(kymax)

c
c     check if k <= kmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

c      leftx = 0
c
c      do 10 i = 1,nx+kx-1
c         if (xknot(i) .gt. xknot(i+1)) then
c             write(6,*) 'subroutine dbs2vl:'
c             write(6,*) 'xknot(i) <= xknot(i+1) required.'
c             write(6,*) i, xknot(i), xknot(i+1)
c             write(6,*)
c             write(6,*) xknot
c          stop
c          endif
c         if((xknot(i) .le. x) .and. (x .lt. xknot(i+1))) leftx = i
c10    continue

cJCH  replaced linear search with binary search
      leftx = ibsearch (x, xknot, nx+kx)

      if(leftx .eq. 0) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'i with xknot(i) <= x < xknot(i+1) required.'
         write(6,*) 'x = ', x
         write(6,*)
         write(6,*) xknot
         stop
      endif

c      lefty = 0

c      do 20 i = 1,ny+ky-1
c         if (yknot(i) .gt. yknot(i+1)) then
c             write(6,*) 'subroutine dbs2vl:'
c             write(6,*) 'yknot(i) <= yknot(i+1) required.'
c             write(6,*) i, yknot(i), yknot(i+1)
c          stop
c          endif
c         if((yknot(i) .le. y) .and. (y .lt. yknot(i+1))) lefty = i
c 20   continue

cJCH  replaced linear search with binary search
      lefty = ibsearch (y, yknot, ny+ky)

      if(lefty .eq. 0) then
         write(6,*) 'subroutine dbs2vl:'
         write(6,*) 'i with yknot(i) <= y < yknot(i+1) required.'
c         write(6,*) 'yknot(i)   = ', yknot(i)
         write(6,*) '  y        = ', y
c         write(6,*) 'yknot(i+1) = ', yknot(i+1)
         stop
      endif
      
      do 30 iky = 1,ky
         work(iky) = dbsdca(0,x,kx,xknot,nx,bcoef(1,lefty-ky+iky),leftx)
 30   continue

      dbs2vl = dbsval(y,ky,yknot(lefty-ky+1),ky,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,bcoef)

c
c  Evaluates the derivative of a two-dimensional tensor-product spline,
c  given its tensor-product B-spline representation.
c
c   iderx  - order of the derivative in the x-direction.  (input)
c   idery  - order of the derivative in the y-direction.  (input)
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   bcoef  - array of length nx*ny containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny.
c   dbs2dr  - value of the (iderx,idery) derivative of the spline at
c            (x,y).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)

      dimension xknot(nx+kx),yknot(ny+ky),bcoef(nx,ny)
      dimension work(kymax)

c
c     check if k <= kmax
c
      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

      nintx = 0

      do 10 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs2dr:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            stop
         endif
         if((xknot(i) .le. x) .and. (x .lt. xknot(i+1))) nintx = i
 10   continue

      if(nintx .eq. 0) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'i with xknot(i) <= x < xknot(i+1) required.'
         write(6,*) 'x = ', x
         stop
      endif

      ninty = 0

      do 20 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
             write(6,*) 'subroutine dbs2dr:'
             write(6,*) 'yknot(i) <= yknot(i+1) required.'
             write(6,*) i, yknot(i), yknot(i+1)
          stop
          endif
         if((yknot(i) .le. y) .and. (y .lt. yknot(i+1))) ninty = i
 20   continue

      if(ninty .eq. 0) then
         write(6,*) 'subroutine dbs2dr:'
         write(6,*) 'i with yknot(i) <= y < yknot(i+1) required.'
         write(6,*) 'y = ', y
         stop
      endif

      do 30 iy = 1, ky
         work(iy) =
     .        dbsdca(iderx,x,kx,xknot,nx,bcoef(1,ninty-ky+iy),nintx)
 30   continue

      dbs2dr = dbsder(idery,y,ky,yknot(ninty-ky+1),ky,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs2gd(iderx,idery,nxvec,xvec,nyvec,yvec,kx,
     .     ky,xknot,yknot,nx,ny,bcoef,val,ldvalue)

c
c  Evaluates the derivative of a two-dimensional tensor-product spline,
c  given its tensor-product B-spline representation on a grid.
c
c   iderx   - order of the derivative in the x-direction.  (input)
c   idery   - order of the derivative in the y-direction.  (input)
c   nxvec   - number of grid points in the x-direction.  (input)
c   xvec    - array of length nx containing the x-coordinates at
c             which the spline is to be evaluated.  (input)
c             the points in xvec should be strictly increasing.
c   nyvec   - number of grid points in the y-direction.  (input)
c   yvec    - array of length ny containing the y-coordinates at
c             which the spline is to be evaluated.  (input)
c             the points in yvec should be strictly increasing.
c   kx      - order of the spline in the x-direction.  (input)
c   ky      - order of the spline in the y-direction.  (input)
c   xknot   - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c   yknot   - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c   nx      - number of B-spline coefficients in the x-direction.
c             (input)
c   ny      - number of B-spline coefficients in the y-direction.
c             (input)
c   bcoef   - array of length nx*ny containing the
c             tensor-product B-spline coefficients.  (input)
c             bscoef is treated internally as a matrix of size nx
c             by ny.
c   val     - value of the (iderx,idery) derivative of the spline on
c             the nx by ny grid.  (output)
c             value(i,j) contains the derivative of the spline at the
c             point (xvec(i),yvec(j)).
c   ldf     - leading dimension of value exactly as specified in the
c             dimension statement of the calling program.  (input)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8)
      parameter(nxmax=200,nymax=200)

      dimension xvec(nxvec),xknot(nx+kx)
      dimension yvec(nyvec),yknot(ny+ky)
      dimension bcoef(nx,ny),val(ldvalue,*)

      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),biaty(nymax,kymax)
      dimension leftx(nxmax),lefty(nymax),term(nxmax)
      dimension work(kymax)

      logical same,next

c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs2gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 10 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 10   continue

      do 20 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 20   continue

      do 30 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 30   continue

c
c     check if ky <= kymax
c

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs2gd:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      lefty(1) = 0

      call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

      do 40 iy = 2,nyvec
         lefty(iy) = lefty(iy-1)
         same = (yknot(lefty(iy)) .le. yvec(iy))
     .        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
         if(.not. same ) then
            lefty(iy) = lefty(iy) + 1
            next      = (yknot(lefty(iy)) .le. yvec(iy))
     .           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
            if (.not. next)
     .           call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
         endif
 40   continue

      do 50 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'yknot(i) <= yknot(i+1) required.'
            write(6,*) i, yknot(i), yknot(i+1)
            write(6,*)
            write(6,*) yknot
            stop
         endif
 50   continue

      do 60 i = 1,nyvec
         if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
            write(6,*) 'subroutine dbs2gd:'
            write(6,*) 'iy with yknot(iy) <= y < yknot(iy+1) required.'
            write(6,*) 'y = ', yvec(i)
            stop
         endif
 60   continue

      if ((iderx .eq. 0) .and. (idery .eq. 0)) then

         do 70 ix = 1,nxvec
            biatx(ix,1) = 1.d0
 70      continue

         do 80 ik = 1,kx-1
            do 90 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 90         continue

            do 100 il = 1,ik
               do 110 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 110           continue
 100        continue

            do 120 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 120        continue
 80      continue

         do 130 iy = 1,nyvec
            biaty(iy,1) = 1.d0
 130      continue

         do 140 ik = 1,ky-1
            do 150 iy = 1,nyvec
               dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
               dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
               save1(iy) = 0.d0
 150         continue

            do 160 il = 1,ik
               do 170 iy = 1,nyvec
                  term(iy)     = biaty(iy,il)
     .                 / (dr(iy,il) + dl(iy,ik+1-il))
                  biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                  save1(iy)    = dl(iy,ik+1-il) * term(iy)
 170           continue
 160        continue

            do 180 iy = 1,nyvec
               biaty(iy,ik+1) = save1(iy)
 180        continue
 140     continue

         do 190 iy = 1,nyvec
            do 200 ix = 1,nxvec
               val(ix,iy) = 0.0d0
 200        continue
 190     continue

         do 210 iky = 1,ky
            do 220 ikx = 1,kx
               do 230 iy = 1,nyvec
                  do 240 ix = 1,nxvec
                     val(ix,iy) = val(ix,iy)
     .                    + biatx(ix,ikx) * biaty(iy,iky)
     .                    * bcoef(leftx(ix)-kx+ikx,lefty(iy)-ky+iky)
 240              continue
 230           continue
 220        continue
 210     continue


      elseif (((iderx .ge. 1) .or. (idery .ge. 1))
     .     .and. ( (iderx .lt. kx) .and. (idery .lt. ky))) then

         do 250 iy = 1,nyvec
            do 260 ix = 1,nxvec
               do 270 iky = 1, ky
                  work(iky) = dbsdca(iderx,xvec(ix),kx,xknot,nx,
     .                 bcoef(1,lefty(iy)-ky+iky),leftx(ix))
 270           continue
               val(ix,iy) = dbsder(idery,yvec(iy),ky,
     .              yknot(lefty(iy)-ky+1),ky,work)
 260        continue
 250     continue

c               val(ix,iy) = dbs2dr(iderx,idery,xvec(ix),yvec(iy),
c     .              kx,ky,xknot,yknot,nx,ny,bcoef)

      else

         do 280 iy = 1,nyvec
            do 290 ix = 1,nxvec
               val(ix,iy) = 0.0d0
 290        continue
 280     continue

      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine dbs3in(nx,xvec,ny,yvec,nz,zvec,xyzdata,ldf,mdf,
     .     kx,ky,kz,xknot,yknot,zknot,bcoef)

c
c  Computes a three-dimensional tensor-product spline interpolant,
c  returning the tensor-product B-spline coefficients.
c
c   nx      - number of data points in the x-direction.  (input)
c   xvec    - array of length nxdata containing the data points in
c             the x-direction.  (input)
c             xdata must be increasing.
c   ny      - number of data points in the y-direction.  (input)
c   yvec    - array of length nydata containing the data points in
c             the y-direction.  (input)
c             ydata must be increasing.
c   nz      - number of data points in the z-direction.  (input)
c   zvec    - array of length nzdata containing the data points in
c             the z-direction.  (input)
c             zdata must be increasing.
c   xyzdata - array of size nx by ny by nz containing the
c             values to be interpolated.  (input)
c             xyzdata(i,j,k) contains the value at
c             (xvec(i),yvec(j),zvec(k)).
c   ldf     - leading dimension of fdata exactly as specified in the
c             dimension statement of the calling program.  (input)
c   mdf     - middle dimension of fdata exactly as specified in the
c             dimension statement of the calling program.  (input)
c   kx      - order of the spline in the x-direction.  (input)
c             kxord must be less than or equal to nxdata.
c   ky      - order of the spline in the y-direction.  (input)
c             kyord must be less than or equal to nydata.
c   kz      - order of the spline in the z-direction.  (input)
c             kzord must be less than or equal to nzdata.
c   xknot   - array of length nx+kx containing the knot
c             sequence in the x-direction.  (input)
c             xknot must be nondecreasing.
c   yknot   - array of length ny+ky containing the knot
c             sequence in the y-direction.  (input)
c             yknot must be nondecreasing.
c   zknot   - array of length nz+kz containing the knot
c             sequence in the z-direction.  (input)
c             zknot must be nondecreasing.
c   bcoef   - array of length nx*ny*nz containing the
c             tensor-product B-spline coefficients.  (output)
c             bscoef is treated internally as a matrix of size nx
c             by ny by nz.
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(nxmax=130,kxmax=8,nymax=130,kymax=8,nzmax=130,kzmax=8)

c
c     dimensions should be
c              work1(nx,ny,nz)
c              work2(nz)
c              work3((2*kz-1)*nz)
c

      dimension work2(nzmax),work3((2*kzmax-1)*nzmax)
      dimension work1(nxmax,nymax,nzmax)

      dimension bcoef(nx,ny,nz),xvec(nx),yvec(ny),zvec(nz)
      dimension xknot(nx+kx),yknot(ny+ky),zknot(nz+kz)
      dimension xyzdata(ldf,mdf,nz)

      if ((kx .gt. kxmax) .or. (nx .gt. nxmax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'kx > kxmax or nx > nxmax'
         write(6,*) 'kx = ',kx,'  kxmax = ',kxmax
         write(6,*) 'nx = ',nx,'  nxmax = ',nxmax
         stop
      endif

      if ((ky .gt. kymax) .or. (ny .gt. nymax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'ky > kymax or ny > nymax'
         write(6,*) 'ky = ',ky,'  kymax = ',kymax
         write(6,*) 'ny = ',ny,'  nymax = ',nymax
         stop
      endif

      if ((kz .gt. kzmax) .or. (nz .gt. nzmax)) then
         write(6,*) 'subroutine dbs3in: error'
         write(6,*) 'kz > kzmax or nz > nzmax'
         write(6,*) 'kz = ',kz,'  kymax = ',kzmax
         write(6,*) 'nz = ',nz,'  nymax = ',nzmax
         stop
      endif

      call spli3d(zvec,ldf,mdf,xyzdata,zknot,nz,kz,nx,ny,work2,work3,
     .   work1,nxmax,nymax,nzmax)

      do 10 iz = 1,nz
         call dbs2in(nx,xvec,ny,yvec,work1(1,1,iz),nxmax,kx,ky,xknot,
     .        yknot,bcoef(1,1,iz))
 10   continue

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine spli3d(xyzvec,ldf,mdf,xyzdata,xyzknot,n,k,m,l,work2,
     .		work3,bcoef,nxmax,nymax,nzmax)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xyzvec(n),xyzknot(n+k),xyzdata(ldf,mdf,*)
      dimension bcoef(nxmax,nymax,nzmax), work2(n),work3((2*k-1)*n)

      np1   = n + 1
      km1   = k - 1
      kpkm2 = 2 * km1
      left  = k
      lenq  = n * (k + km1)

      do 10 i = 1,lenq
         work3(i) = 0.d0
 10   continue

      do 20 i = 1,n
         xyzveci = xyzvec(i)
         ilp1mx  = min0(i+k,np1)
         left    = max0(left,i)
         if (xyzveci .lt. xyzknot(left)) go to 998
 30      if (xyzveci .lt. xyzknot(left+1)) go to 40
         left = left + 1
         if (left .lt. ilp1mx) go to 30
         left = left - 1
         if (xyzveci .gt. xyzknot(left+1)) go to 998
 40      call bsplvb(xyzknot,n+k,k,1,xyzveci,left,work2)
         jj = i - left + 1 + (left - k) * (k + km1)
         do 50 j = 1,k
            jj    = jj + kpkm2
            work3(jj) = work2(j)
 50      continue
 20   continue

      call banfac(work3,k+km1,n,km1,km1,iflag)

      go to (60,999), iflag

 60   do 70 j = 1,l
         do 80 i = 1,m
            do 90 in = 1,n
               work2(in) = xyzdata(i,j,in)
 90         continue

            call banslv(work3,k+km1,n,km1,km1,work2)

            do 100 in = 1,n
               bcoef(i,j,in) = work2(in)
 100        continue
 80      continue
 70   continue

      return

 998  write(6,*) 'subroutine db3in:'
      write(6,*) 'i with knot(i) <= x/y/z < knot(i+1) required.'
      write(6,*) 'knot(1)   = ', xyzknot(1)
      write(6,*) 'knot(n+k) = ', xyzknot(n+k)
      write(6,*) '    x/y/z = ', xyzveci

      stop

 999  write(6,*) 'subroutine dbs3in: error'
      write(6,*) 'no solution of linear equation system !!!'

      stop

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs3vl(x,y,z,kx,ky,kz,xknot,yknot,zknot,nx,ny,nz,bcoef)

c
c  Evaluates a three-dimensional tensor-product spline, given its
c  tensor-product B-spline representation.
c
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   z      - z-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   dbs3vl - value of the spline at (x,y,z).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)

c
c     dimension should be
c            dimension work(kz)
c

      dimension work(kzmax)

      dimension xknot(nx+kx),yknot(ny+ky),zknot(nz+kz)
      dimension bcoef(nx,ny,nz)

c     check if k <= kmax

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c


      nintz = 0

      do 10 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i + 1)) then
             write(6,*) 'subroutine dbs3vl:'
             write(6,*) 'zknot(i) <= zknot(i+1) required.'
             write(6,*) i, zknot(i), zknot(i+1)
          stop
          endif
         if((zknot(i) .le. z) .and. (z .lt. zknot(i + 1))) nintz = i
 10   continue

      if(nintz .eq. 0) then
         write(6,*) 'subroutine dbs3vl:'
         write(6,*) 'i with zknot(i) <= z < zknot(i+1) required.'
         write(6,*) 'zknot(i)   = ', zknot(i)
         write(6,*) '  z        = ', z
         write(6,*) 'zknot(i+1) = ', zknot(i+1)
         stop
      endif

      do 40 jz = 1,kz
         work(jz) = dbs2vl(x,y,kx,ky,xknot,yknot,nx,ny,
     .        bcoef(1,1,nintz-kz+jz))
 40   continue

      dbs3vl = dbsval(z,kz,zknot(nintz-kz+1),kz,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      function dbs3dr(iderx,idery,iderz,x,y,z,kx,ky,kz,
     .     xknot,yknot,zknot,nx,ny,nz,bcoef)

c
c  Evaluates the derivative of a three-dimensional tensor-product spline,
c  given its tensor-product B-spline representation.
c
c   iderx  - order of the x-derivative.  (input)
c   idery  - order of the y-derivative.  (input)
c   iderz  - order of the z-derivative.  (input)
c   x      - x-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   y      - y-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   z      - z-coordinate of the point at which the spline is to be
c            evaluated.  (input)
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   dbs3dr - value of the (iderx,idery,iderz) derivative of the
c            spline at (x,y,z).  (output)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)

      dimension xknot(nx+kx)
      dimension yknot(ny+ky)
      dimension zknot(nz+kz)
      dimension bcoef(nx,ny,nz),work(kzmax)

c
c     check if k <= kmax
c
      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

c
c     check if knot(i) <= knot(i+1) and calculation of i so that
c     knot(i) <= x < knot(i+1)
c

      nintz = 0

      do 10 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i + 1)) then
             write(6,*) 'subroutine dbs3vl:'
             write(6,*) 'zknot(i) <= zknot(i+1) required.'
             write(6,*) i, zknot(i), zknot(i+1)
          stop
          endif
         if((zknot(i) .le. z) .and. (z .lt. zknot(i + 1))) nintz = i
 10   continue

      if(nintz .eq. 0) then
         write(6,*) 'subroutine dbs3dr:'
         write(6,*) 'i with zknot(i) <= z < zknot(i+1) required.'
         write(6,*) 'zknot(i)   = ', zknot(i)
         write(6,*) '  z        = ', z
         write(6,*) 'zknot(i+1) = ', zknot(i+1)
         stop
      endif

      do 20 jz = 1,kz
         work(jz) = dbs2dr(iderx,idery,x,y,kx,ky,xknot,yknot,nx,ny,
     .        bcoef(1,1,nintz-kz+jz))
 20   continue

      dbs3dr = dbsder(iderz,z,kz,zknot(nintz-kz+1),kz,work)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine dbs3gd(iderx,idery,iderz,nxvec,xvec,nyvec,yvec,nzvec,
     .     zvec,kx,ky,kz,xknot,yknot,zknot,nx,ny,
     .     nz,bcoef,value,ldvalue,mdvalue)

c
c  Evaluates the derivative of a three-dimensional tensor-product spline,
c  given its tensor-product B-spline representation on a grid.
c
c   iderx  - order of the x-derivative.  (input)
c   idery  - order of the y-derivative.  (input)
c   iderz  - order of the z-derivative.  (input)
c   nx     - number of grid points in the x-direction.  (input)
c   xvec   - array of length nx containing the x-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in xvec should be strictly increasing.
c   ny     - number of grid points in the y-direction.  (input)
c   yvec   - array of length ny containing the y-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in yvec should be strictly increasing.
c   nz     - number of grid points in the z-direction.  (input)
c   zvec   - array of length nz containing the z-coordinates at
c            which the spline is to be evaluated.  (input)
c            the points in yvec should be strictly increasing.
c   kx     - order of the spline in the x-direction.  (input)
c   ky     - order of the spline in the y-direction.  (input)
c   kz     - order of the spline in the z-direction.  (input)
c   xknot  - array of length nx+kx containing the knot
c            sequence in the x-direction.  (input)
c            xknot must be nondecreasing.
c   yknot  - array of length ny+ky containing the knot
c            sequence in the y-direction.  (input)
c            yknot must be nondecreasing.
c   zknot  - array of length nz+kz containing the knot
c            sequence in the z-direction.  (input)
c            zknot must be nondecreasing.
c   nx     - number of B-spline coefficients in the x-direction.
c            (input)
c   ny     - number of B-spline coefficients in the y-direction.
c            (input)
c   nz     - number of B-spline coefficients in the z-direction.
c            (input)
c   bcoef  - array of length nx*ny*nz containing the
c            tensor-product B-spline coefficients.  (input)
c            bscoef is treated internally as a matrix of size nx
c            by ny by nz.
c   val    - array of size nx by ny by nz containing the values of
c            the (iderx,idery,iderz) derivative of the spline on the
c            nx by ny by nz grid.  (output)
c            value(i,j,k) contains the derivative of the spline at
c            the point (xvec(i), yvec(j), zvec(k)).
c   ldf    - leading dimension of value exactly as specified in the
c            dimension statement of the calling program.  (input)
c   mdf    - middle dimension of value exactly as specified in the
c            dimension statement of the calling program.  (input)
c

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kxmax=8,kymax=8,kzmax=8)
      parameter(nxmax=200,nymax=200,nzmax=200)

      dimension xvec(nxvec),xknot(nx+kx)
      dimension yvec(nyvec),yknot(ny+ky)
      dimension zvec(nzvec),zknot(nz+kz)
      dimension bcoef(nx,ny,nz)
      dimension value(ldvalue,mdvalue,*)

      dimension dl(nxmax,kxmax),dr(nxmax,kxmax),save1(nxmax)
      dimension biatx(nxmax,kxmax),biaty(nymax,kymax)
      dimension biatz(nzmax,kzmax)
      dimension leftx(nxmax),lefty(nymax),leftz(nzmax)
      dimension term(nxmax)

      logical same,next

c
c     check if kx <= kxmax
c

      if(kx .gt. kxmax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'kx <= kxmax required.'
         stop
      endif

      do 10 i = 1,nx+kx-1
         if (xknot(i) .gt. xknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'xknot(i) <= xknot(i+1) required.'
            write(6,*) i, xknot(i), xknot(i+1)
            write(6,*)
            write(6,*) xknot
            stop
         endif
 10   continue

      do 20 i = 1,nxvec
         if ((xvec(i).lt.xknot(1)).or.(xvec(i).gt.xknot(nx+kx))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'ix with xknot(ix) <= x < xknot(ix+1) required.'
            write(6,*) 'x = ', xvec(i)
            stop
         endif
 20   continue

      leftx(1) = 0

      call huntn(xknot,nx+kx,kx,xvec(1),leftx(1))

      do 30 ix = 2,nxvec
         leftx(ix) = leftx(ix-1)
         same = (xknot(leftx(ix)) .le. xvec(ix))
     .        .and. (xvec(ix) .le. xknot(leftx(ix)+1))
         if(.not. same ) then
            leftx(ix) = leftx(ix) + 1
            next      = (xknot(leftx(ix)) .le. xvec(ix))
     .           .and. (xvec(ix) .le. xknot(leftx(ix)+1))
            if (.not. next)
     .           call huntn(xknot,nx+kx,kx,xvec(ix),leftx(ix))
         endif
 30   continue

c
c     check if ky <= kymax
c

      if(ky .gt. kymax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'ky <= kymax required.'
         stop
      endif

      do 40 i = 1,ny+ky-1
         if (yknot(i) .gt. yknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'yknot(i) <= yknot(i+1) required.'
            write(6,*) i, yknot(i), yknot(i+1)
            write(6,*)
            write(6,*) yknot
            stop
         endif
 40   continue

      do 50 i = 1,nyvec
         if ((yvec(i).lt.yknot(1)).or.(yvec(i).gt.yknot(ny+ky))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'iy with yknot(iy) <= y < yknot(iy+1) required.'
            write(6,*) 'y = ', yvec(i)
            stop
         endif
 50   continue

      lefty(1) = 0

      call huntn(yknot,ny+ky,ky,yvec(1),lefty(1))

      do 60 iy = 2,nyvec
         lefty(iy) = lefty(iy-1)
         same = (yknot(lefty(iy)) .le. yvec(iy))
     .        .and. (yvec(iy) .le. yknot(lefty(iy)+1))
         if(.not. same ) then
            lefty(iy) = lefty(iy) + 1
            next      = (yknot(lefty(iy)) .le. yvec(iy))
     .           .and. (yvec(iy) .le. yknot(lefty(iy)+1))
            if (.not. next)
     .           call huntn(yknot,ny+ky,ky,yvec(iy),lefty(iy))
         endif
 60   continue

c
c     check if kz <= kzmax
c

      if(kz .gt. kzmax) then
         write(6,*) 'subroutine dbs3gd:'
         write(6,*) 'kz <= kzmax required.'
         stop
      endif

      do 70 i = 1,nz+kz-1
         if (zknot(i) .gt. zknot(i+1)) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'zknot(i) <= zknot(i+1) required.'
            write(6,*) i, zknot(i), zknot(i+1)
            write(6,*)
            write(6,*) zknot
            stop
         endif
 70   continue

      do 80 i = 1,nzvec
         if ((zvec(i).lt.zknot(1)).or.(zvec(i).gt.zknot(nz+kz))) then
            write(6,*) 'subroutine dbs3gd:'
            write(6,*) 'iz with zknot(iz) <= z < zknot(iz+1) required.'
            write(6,*) 'z = ', zvec(i)
            stop
         endif
 80   continue

      leftz(1) = 0

      call huntn(zknot,nz+kz,kz,zvec(1),leftz(1))

      do 90 iz = 2,nzvec
         leftz(iz) = leftz(iz-1)
         same = (zknot(leftz(iz)) .le. zvec(iz))
     .        .and. (zvec(iz) .le. zknot(leftz(iz)+1))
         if(.not. same ) then
            leftz(iz) = leftz(iz) + 1
            next      = (zknot(leftz(iz)) .le. zvec(iz))
     .           .and. (zvec(iz) .le. zknot(leftz(iz)+1))
            if (.not. next)
     .           call huntn(zknot,nz+kz,kz,zvec(iz),leftz(iz))
         endif
 90   continue

      if ((iderx .eq. 0) .and. (idery .eq. 0) .and. (iderz .eq.0)) then

         do 100 ix = 1,nxvec
            biatx(ix,1) = 1.d0
 100      continue

         do 110 ik = 1,kx-1
            do 120 ix = 1,nxvec
               dr(ix,ik) = xknot(leftx(ix)+ik) - xvec(ix)
               dl(ix,ik) = xvec(ix) - xknot(leftx(ix)+1-ik)
               save1(ix) = 0.d0
 120         continue

            do 130 il = 1,ik
               do 140 ix = 1,nxvec
                  term(ix)     = biatx(ix,il)
     .                 / (dr(ix,il) + dl(ix,ik+1-il))
                  biatx(ix,il) = save1(ix) + dr(ix,il) * term(ix)
                  save1(ix)    = dl(ix,ik+1-il) * term(ix)
 140           continue
 130        continue

            do 150 ix = 1,nxvec
               biatx(ix,ik+1) = save1(ix)
 150        continue
 110     continue

         do 160 iy = 1,nyvec
            biaty(iy,1) = 1.d0
 160     continue

         do 170 ik = 1,ky-1
            do 180 iy = 1,nyvec
               dr(iy,ik) = yknot(lefty(iy)+ik) - yvec(iy)
               dl(iy,ik) = yvec(iy) - yknot(lefty(iy)+1-ik)
               save1(iy) = 0.d0
 180        continue

            do 190 il = 1,ik
               do 200 iy = 1,nyvec
                  term(iy)     = biaty(iy,il)
     .                 / (dr(iy,il) + dl(iy,ik+1-il))
                  biaty(iy,il) = save1(iy) + dr(iy,il) * term(iy)
                  save1(iy)    = dl(iy,ik+1-il) * term(iy)
 200           continue
 190        continue

            do 210 iy = 1,nyvec
               biaty(iy,ik+1) = save1(iy)
 210        continue
 170     continue

         do 220 iz = 1,nzvec
            biatz(iz,1) = 1.d0
 220     continue

         do 230 ik = 1,kz-1
            do 240 iz = 1,nzvec
               dr(iz,ik) = zknot(leftz(iz)+ik) - zvec(iz)
               dl(iz,ik) = zvec(iz) - zknot(leftz(iz)+1-ik)
               save1(iz) = 0.d0
 240        continue

            do 250 il = 1,ik
               do 260 iz = 1,nzvec
                  term(iz)     = biatz(iz,il)
     .                 / (dr(iz,il) + dl(iz,ik+1-il))
                  biatz(iz,il) = save1(iz) + dr(iz,il) * term(iz)
                  save1(iz)    = dl(iz,ik+1-il) * term(iz)
 260           continue
 250        continue

            do 270 iz = 1,nzvec
               biatz(iz,ik+1) = save1(iz)
 270        continue
 230     continue


         do 280 iz = 1,nzvec
            do 290 iy = 1,nyvec
               do 300 ix = 1,nxvec
                  value(ix,iy,iz) = 0.0d0
 300           continue
 290        continue
 280     continue

         do 310 ikz = 1,kz
            do 320 iky = 1,ky
               do 330 ikx = 1,kx
                  do 340 iz = 1,nzvec
                     do 350 iy = 1,nyvec
                        do 360 ix = 1,nxvec
                           value(ix,iy,iz) = value(ix,iy,iz)
     .                          + biatx(ix,ikx) * biaty(iy,iky)
     .                          * biatz(iz,ikz)
     .                          * bcoef(leftx(ix)-kx+ikx,
     .                          lefty(iy)-ky+iky,leftz(iz)-kz+ikz)
 360                    continue
 350                 continue
 340              continue
 330           continue
 320        continue
 310     continue

      else

         do 370 iz = 1,nzvec
            do 380 iy = 1,nyvec
               do 390 ix = 1,nxvec
                  value(ix,iy,iz) = dbs3dr(iderx,idery,iderz,xvec(ix),
     .                 yvec(iy),zvec(iz),kx,ky,kz,xknot,yknot,
     .                 zknot,nx,ny,nz,bcoef)
 390           continue
 380         continue
 370      continue

      endif

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

      implicit double precision (a-h,o-z), integer (i-n)

      parameter(kmax=8)

      dimension biatx(jhigh),t(n),dl(kmax),dr(kmax)

      data j/1/

      go to (10,20), index

 10   j = 1

      biatx(1) = 1.d0

      if (j .ge. jhigh) go to 99

 20   jp1 = j + 1

      dr(j) = t(left+j) - x
      dl(j) = x - t(left+1-j)
      saved = 0.d0

      do 30 i = 1, j
         term     = biatx(i) / (dr(i) + dl(jp1-i))
         biatx(i) = saved + dr(i) * term
         saved    = dl(jp1-i) * term
 30   continue

      biatx(jp1) = saved
      j          = jp1

      if (j .lt. jhigh) go to 20

 99   return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow)

      iflag  = 1
      middle = nbandu + 1
      nrowm1 = nrow - 1

      if (nrowm1) 999,900,10

 10   if (nbandl .gt. 0) go to 30

      do 20 i = 1,nrowm1
         if (w(middle,i) .eq. 0.d0) go to 999
 20   continue

      go to 900

 30   if (nbandu .gt. 0) go to 60
      do 40 i = 1,nrowm1
         pivot = w(middle,i)
         if(pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl, nrow - i)
         do 50 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 50      continue
 40   continue

      return

 60   do 70 i = 1,nrowm1
         pivot = w(middle,i)
         if (pivot .eq. 0.d0) go to 999
         jmax = min0(nbandl,nrow - i)
         do 80 j = 1,jmax
            w(middle+j,i) = w(middle+j,i) / pivot
 80      continue

         kmax = min0(nbandu,nrow - i)

         do 90 k = 1,kmax
            ipk    = i + k
            midmk  = middle - k
            factor = w(midmk,ipk)
            do 100 j = 1,jmax
               w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)
     .              * factor
 100         continue
 90       continue
 70   continue

 900  if (w(middle,nrow) .ne. 0.d0) return
 999  iflag = 2

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension w(nroww,nrow),b(nrow)

      middle = nbandu + 1
      if (nrow .eq. 1) goto 99
      nrowm1 = nrow - 1
      if (nbandl .eq. 0) goto 30

      do 10 i = 1, nrowm1
         jmax = min0(nbandl, nrow - i)
         do 20 j = 1, jmax
            b(i+j) = b(i+j) - b(i) * w(middle+j,i)
 20      continue
 10   continue

 30   if (nbandu .gt. 0)  goto 50
      do 40 i = 1, nrow
         b(i) = b(i) / w(1,i)
 40   continue

      return

 50   do 60 i = nrow, 2, -1
         b(i) = b(i)/w(middle,i)
         jmax = min0(nbandu,i-1)
         do 70 j = 1,jmax
            b(i-j) = b(i-j) - b(i) * w(middle-j,i)
 70      continue
 60   continue

 99   b(1) = b(1) / w(middle,1)

      return

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      subroutine huntn(xx,n,kord,x,jlo)

      implicit double precision (a-h,o-z), integer (i-n)

      dimension xx(n)

c
c     works only for B-Splines (order n)
c

      max  = n - kord
      null = kord

      if (jlo.le.null.or.jlo.gt.max) then
         jlo = null
         jhi = max+1
         goto 30
      endif

      inc = 1

      if (x .ge. xx(jlo)) then
 10      jhi = jlo + inc
         if (jhi .gt. max) then
            jhi = max + 1
         else if (x .ge. xx(jhi)) then
            jlo = jhi
            inc = inc + inc
            goto 10
         endif
      else
         jhi = jlo
 20      jlo = jhi - inc
         if (jlo .le. null) then
            jlo = null
         else if (x .lt. xx(jlo)) then
            jhi = jlo
            inc = inc + inc
            goto 20
         endif
      endif

 30   if (jhi-jlo.eq.1) return

      jm = (jhi + jlo) / 2
      if (x .gt. xx(jm)) then
         jlo = jm
      else
         jhi = jm
      endif

      goto 30

      end


c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
