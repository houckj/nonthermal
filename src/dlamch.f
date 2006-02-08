      double precision function dlamch( cmach )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      character          cmach
*     ..
*
*  purpose
*  =======
*
*  dlamch determines double precision machine parameters.
*
*  arguments
*  =========
*
*  cmach   (input) character*1
*          specifies the value to be returned by dlamch:
*          = 'e' or 'e',   dlamch := eps
*          = 's' or 's ,   dlamch := sfmin
*          = 'b' or 'b',   dlamch := base
*          = 'p' or 'p',   dlamch := eps*base
*          = 'n' or 'n',   dlamch := t
*          = 'r' or 'r',   dlamch := rnd
*          = 'm' or 'm',   dlamch := emin
*          = 'u' or 'u',   dlamch := rmin
*          = 'l' or 'l',   dlamch := emax
*          = 'o' or 'o',   dlamch := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      logical            first, lrnd
      integer            beta, imax, imin, it
      double precision   base, emax, emin, eps, prec, rmach, rmax, rmin,
     $                   rnd, sfmin, small, t
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           dlamc2
*     ..
*     .. save statement ..
      save               first, eps, sfmin, base, t, rnd, emin, rmin,
     $                   emax, rmax, prec
*     ..
*     .. data statements ..
      data               first / .true. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         call dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         if( lrnd ) then
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         else
            rnd = zero
            eps = base**( 1-it )
         end if
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         if( small.ge.sfmin ) then
*
*           use small plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         end if
      end if
*
      if( lsame( cmach, 'e' ) ) then
         rmach = eps
      else if( lsame( cmach, 's' ) ) then
         rmach = sfmin
      else if( lsame( cmach, 'b' ) ) then
         rmach = base
      else if( lsame( cmach, 'p' ) ) then
         rmach = prec
      else if( lsame( cmach, 'n' ) ) then
         rmach = t
      else if( lsame( cmach, 'r' ) ) then
         rmach = rnd
      else if( lsame( cmach, 'm' ) ) then
         rmach = emin
      else if( lsame( cmach, 'u' ) ) then
         rmach = rmin
      else if( lsame( cmach, 'l' ) ) then
         rmach = emax
      else if( lsame( cmach, 'o' ) ) then
         rmach = rmax
      end if
*
      dlamch = rmach
      return
*
*     end of dlamch
*
      end
*
************************************************************************
*
      subroutine dlamc1( beta, t, rnd, ieee1 )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            ieee1, rnd
      integer            beta, t
*     ..
*
*  purpose
*  =======
*
*  dlamc1 determines the machine parameters given by beta, t, rnd, and
*  ieee1.
*
*  arguments
*  =========
*
*  beta    (output) integer
*          the base of the machine.
*
*  t       (output) integer
*          the number of ( beta ) digits in the mantissa.
*
*  rnd     (output) logical
*          specifies whether proper rounding  ( rnd = .true. )  or
*          chopping  ( rnd = .false. )  occurs in addition. this may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  ieee1   (output) logical
*          specifies whether rounding appears to be done in the ieee
*          'round to nearest' style.
*
*  further details
*  ===============
*
*  the routine is based on the routine  envron  by malcolm and
*  incorporates suggestions by gentleman and marovich. see
*
*     malcolm m. a. (1972) algorithms to reveal properties of
*        floating-point arithmetic. comms. of the acm, 15, 949-951.
*
*     gentleman w. m. and marovich s. b. (1974) more on algorithms
*        that reveal properties of floating point arithmetic units.
*        comms. of the acm, 17, 276-277.
*
* =====================================================================
*
*     .. local scalars ..
      logical            first, lieee1, lrnd
      integer            lbeta, lt
      double precision   a, b, c, f, one, qtr, savec, t1, t2
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. save statement ..
      save               first, lieee1, lbeta, lrnd, lt
*     ..
*     .. data statements ..
      data               first / .true. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         one = 1
*
*        lbeta,  lieee1,  lt and  lrnd  are the  local values  of  beta,
*        ieee1, t and rnd.
*
*        throughout this routine  we use the function  dlamc3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         a = 1
         c = 1
*
*+       while( c.eq.one )loop
   10    continue
         if( c.eq.one ) then
            a = 2*a
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            go to 10
         end if
*+       end while
*
*        now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         b = 1
         c = dlamc3( a, b )
*
*+       while( c.eq.a )loop
   20    continue
         if( c.eq.a ) then
            b = 2*b
            c = dlamc3( a, b )
            go to 20
         end if
*+       end while
*
*        now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         qtr = one / 4
         savec = c
         c = dlamc3( c, -a )
         lbeta = c + qtr
*
*        now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         b = lbeta
         f = dlamc3( b / 2, -b / 100 )
         c = dlamc3( f, a )
         if( c.eq.a ) then
            lrnd = .true.
         else
            lrnd = .false.
         end if
         f = dlamc3( b / 2, b / 100 )
         c = dlamc3( f, a )
         if( ( lrnd ) .and. ( c.eq.a ) )
     $      lrnd = .false.
*
*        try and decide whether rounding is done in the  ieee  'round to
*        nearest' style. b/2 is half a unit in the last place of the two
*        numbers a and savec. furthermore, a is even, i.e. has last  bit
*        zero, and savec is odd. thus adding b/2 to a should not  change
*        a, but adding b/2 to savec should change savec.
*
         t1 = dlamc3( b / 2, a )
         t2 = dlamc3( b / 2, savec )
         lieee1 = ( t1.eq.a ) .and. ( t2.gt.savec ) .and. lrnd
*
*        now find  the  mantissa, t.  it should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  so we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         lt = 0
         a = 1
         c = 1
*
*+       while( c.eq.one )loop
   30    continue
         if( c.eq.one ) then
            lt = lt + 1
            a = a*lbeta
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            go to 30
         end if
*+       end while
*
      end if
*
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      return
*
*     end of dlamc1
*
      end
*
************************************************************************
*
      subroutine dlamc2( beta, t, rnd, eps, emin, rmin, emax, rmax )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            rnd
      integer            beta, emax, emin, t
      double precision   eps, rmax, rmin
*     ..
*
*  purpose
*  =======
*
*  dlamc2 determines the machine parameters specified in its argument
*  list.
*
*  arguments
*  =========
*
*  beta    (output) integer
*          the base of the machine.
*
*  t       (output) integer
*          the number of ( beta ) digits in the mantissa.
*
*  rnd     (output) logical
*          specifies whether proper rounding  ( rnd = .true. )  or
*          chopping  ( rnd = .false. )  occurs in addition. this may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  eps     (output) double precision
*          the smallest positive number such that
*
*             fl( 1.0 - eps ) .lt. 1.0,
*
*          where fl denotes the computed value.
*
*  emin    (output) integer
*          the minimum exponent before (gradual) underflow occurs.
*
*  rmin    (output) double precision
*          the smallest normalized number for the machine, given by
*          base**( emin - 1 ), where  base  is the floating point value
*          of beta.
*
*  emax    (output) integer
*          the maximum exponent before overflow occurs.
*
*  rmax    (output) double precision
*          the largest positive number for the machine, given by
*          base**emax * ( 1 - eps ), where  base  is the floating point
*          value of beta.
*
*  further details
*  ===============
*
*  the computation of  eps  is based on a routine paranoia by
*  w. kahan of the university of california at berkeley.
*
* =====================================================================
*
*     .. local scalars ..
      logical            first, ieee, iwarn, lieee1, lrnd
      integer            gnmin, gpmin, i, lbeta, lemax, lemin, lt,
     $                   ngnmin, ngpmin
      double precision   a, b, c, half, leps, lrmax, lrmin, one, rbase,
     $                   sixth, small, third, two, zero
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. external subroutines ..
      external           dlamc1, dlamc4, dlamc5
*     ..
*     .. intrinsic functions ..
      intrinsic          abs, max, min
*     ..
*     .. save statement ..
      save               first, iwarn, lbeta, lemax, lemin, leps, lrmax,
     $                   lrmin, lt
*     ..
*     .. data statements ..
      data               first / .true. / , iwarn / .false. /
*     ..
*     .. executable statements ..
*
      if( first ) then
         first = .false.
         zero = 0
         one = 1
         two = 2
*
*        lbeta, lt, lrnd, leps, lemin and lrmin  are the local values of
*        beta, t, rnd, eps, emin and rmin.
*
*        throughout this routine  we use the function  dlamc3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        dlamc1 returns the parameters  lbeta, lt, lrnd and lieee1.
*
         call dlamc1( lbeta, lt, lrnd, lieee1 )
*
*        start to find eps.
*
         b = lbeta
         a = b**( -lt )
         leps = a
*
*        try some tricks to see whether or not this is the correct  eps.
*
         b = two / 3
         half = one / 2
         sixth = dlamc3( b, -half )
         third = dlamc3( sixth, sixth )
         b = dlamc3( third, -half )
         b = dlamc3( b, sixth )
         b = abs( b )
         if( b.lt.leps )
     $      b = leps
*
         leps = 1
*
*+       while( ( leps.gt.b ).and.( b.gt.zero ) )loop
   10    continue
         if( ( leps.gt.b ) .and. ( b.gt.zero ) ) then
            leps = b
            c = dlamc3( half*leps, ( two**5 )*( leps**2 ) )
            c = dlamc3( half, -c )
            b = dlamc3( half, c )
            c = dlamc3( half, -b )
            b = dlamc3( half, c )
            go to 10
         end if
*+       end while
*
         if( a.lt.leps )
     $      leps = a
*
*        computation of eps complete.
*
*        now find  emin.  let a = + or - 1, and + or - (1 + base**(-3)).
*        keep dividing  a by beta until (gradual) underflow occurs. this
*        is detected when we cannot recover the previous a.
*
         rbase = one / lbeta
         small = one
         do 20 i = 1, 3
            small = dlamc3( small*rbase, zero )
   20    continue
         a = dlamc3( one, small )
         call dlamc4( ngpmin, one, lbeta )
         call dlamc4( ngnmin, -one, lbeta )
         call dlamc4( gpmin, a, lbeta )
         call dlamc4( gnmin, -a, lbeta )
         ieee = .false.
*
         if( ( ngpmin.eq.ngnmin ) .and. ( gpmin.eq.gnmin ) ) then
            if( ngpmin.eq.gpmin ) then
               lemin = ngpmin
*            ( non twos-complement machines, no gradual underflow;
*              e.g.,  vax )
            else if( ( gpmin-ngpmin ).eq.3 ) then
               lemin = ngpmin - 1 + lt
               ieee = .true.
*            ( non twos-complement machines, with gradual underflow;
*              e.g., ieee standard followers )
            else
               lemin = min( ngpmin, gpmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else if( ( ngpmin.eq.gpmin ) .and. ( ngnmin.eq.gnmin ) ) then
            if( abs( ngpmin-ngnmin ).eq.1 ) then
               lemin = max( ngpmin, ngnmin )
*            ( twos-complement machines, no gradual underflow;
*              e.g., cyber 205 )
            else
               lemin = min( ngpmin, ngnmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else if( ( abs( ngpmin-ngnmin ).eq.1 ) .and.
     $            ( gpmin.eq.gnmin ) ) then
            if( ( gpmin-min( ngpmin, ngnmin ) ).eq.3 ) then
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*            ( twos-complement machines with gradual underflow;
*              no known machine )
            else
               lemin = min( ngpmin, ngnmin )
*            ( a guess; no known machine )
               iwarn = .true.
            end if
*
         else
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*         ( a guess; no known machine )
            iwarn = .true.
         end if
***
* comment out this if block if emin is ok
         if( iwarn ) then
            first = .true.
            write( 6, fmt = 9999 )lemin
         end if
***
*
*        assume ieee arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  ieee style,  determined
*        in routine dlamc1. a true ieee machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         ieee = ieee .or. lieee1
*
*        compute  rmin by successive division by  beta. we could compute
*        rmin as base**( emin - 1 ),  but some machines underflow during
*        this computation.
*
         lrmin = 1
         do 30 i = 1, 1 - lemin
            lrmin = dlamc3( lrmin*rbase, zero )
   30    continue
*
*        finally, call dlamc5 to compute emax and rmax.
*
         call dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax )
      end if
*
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
*
      return
*
 9999 format( / / ' warning. the value emin may be incorrect:-',
     $      '  emin = ', i8, /
     $      ' if, after inspection, the value emin looks',
     $      ' acceptable please comment out ',
     $      / ' the if block as marked within the code of routine',
     $      ' dlamc2,', / ' otherwise supply emin explicitly.', / )
*
*     end of dlamc2
*
      end
*
************************************************************************
*
      double precision function dlamc3( a, b )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      double precision   a, b
*     ..
*
*  purpose
*  =======
*
*  dlamc3  is intended to force  a  and  b  to be stored prior to doing
*  the addition of  a  and  b ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  arguments
*  =========
*
*  a, b    (input) double precision
*          the values a and b.
*
* =====================================================================
*
*     .. executable statements ..
*
      dlamc3 = a + b
*
      return
*
*     end of dlamc3
*
      end
*
************************************************************************
*
      subroutine dlamc4( emin, start, base )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      integer            base, emin
      double precision   start
*     ..
*
*  purpose
*  =======
*
*  dlamc4 is a service routine for dlamc2.
*
*  arguments
*  =========
*
*  emin    (output) emin
*          the minimum exponent before (gradual) underflow, computed by
*          setting a = start and dividing by base until the previous a
*          can not be recovered.
*
*  start   (input) double precision
*          the starting point for determining emin.
*
*  base    (input) integer
*          the base of the machine.
*
* =====================================================================
*
*     .. local scalars ..
      integer            i
      double precision   a, b1, b2, c1, c2, d1, d2, one, rbase, zero
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. executable statements ..
*
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = dlamc3( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
*+    while( ( c1.eq.a ).and.( c2.eq.a ).and.
*    $       ( d1.eq.a ).and.( d2.eq.a )      )loop
   10 continue
      if( ( c1.eq.a ) .and. ( c2.eq.a ) .and. ( d1.eq.a ) .and.
     $    ( d2.eq.a ) ) then
         emin = emin - 1
         a = b1
         b1 = dlamc3( a / base, zero )
         c1 = dlamc3( b1*base, zero )
         d1 = zero
         do 20 i = 1, base
            d1 = d1 + b1
   20    continue
         b2 = dlamc3( a*rbase, zero )
         c2 = dlamc3( b2 / rbase, zero )
         d2 = zero
         do 30 i = 1, base
            d2 = d2 + b2
   30    continue
         go to 10
      end if
*+    end while
*
      return
*
*     end of dlamc4
*
      end
*
************************************************************************
*
      subroutine dlamc5( beta, p, emin, ieee, emax, rmax )
*
*  -- lapack auxiliary routine (version 2.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     october 31, 1992
*
*     .. scalar arguments ..
      logical            ieee
      integer            beta, emax, emin, p
      double precision   rmax
*     ..
*
*  purpose
*  =======
*
*  dlamc5 attempts to compute rmax, the largest machine floating-point
*  number, without overflow.  it assumes that emax + abs(emin) sum
*  approximately to a power of 2.  it will fail on machines where this
*  assumption does not hold, for example, the cyber 205 (emin = -28625,
*  emax = 28718).  it will also fail if the value supplied for emin is
*  too large (i.e. too close to zero), probably with overflow.
*
*  arguments
*  =========
*
*  beta    (input) integer
*          the base of floating-point arithmetic.
*
*  p       (input) integer
*          the number of base beta digits in the mantissa of a
*          floating-point value.
*
*  emin    (input) integer
*          the minimum exponent before (gradual) underflow.
*
*  ieee    (input) logical
*          a logical flag specifying whether or not the arithmetic
*          system is thought to comply with the ieee standard.
*
*  emax    (output) integer
*          the largest exponent before overflow
*
*  rmax    (output) double precision
*          the largest machine floating-point number.
*
* =====================================================================
*
*     .. parameters ..
      double precision   zero, one
      parameter          ( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. local scalars ..
      integer            exbits, expsum, i, lexp, nbits, try, uexp
      double precision   oldy, recbas, y, z
*     ..
*     .. external functions ..
      double precision   dlamc3
      external           dlamc3
*     ..
*     .. intrinsic functions ..
      intrinsic          mod
*     ..
*     .. executable statements ..
*
*     first compute lexp and uexp, two powers of 2 that bound
*     abs(emin). we then assume that emax + abs(emin) will sum
*     approximately to the bound that is closest to abs(emin).
*     (emax is the exponent of the required number rmax).
*
      lexp = 1
      exbits = 1
   10 continue
      try = lexp*2
      if( try.le.( -emin ) ) then
         lexp = try
         exbits = exbits + 1
         go to 10
      end if
      if( lexp.eq.-emin ) then
         uexp = lexp
      else
         uexp = try
         exbits = exbits + 1
      end if
*
*     now -lexp is less than or equal to emin, and -uexp is greater
*     than or equal to emin. exbits is the number of bits needed to
*     store the exponent.
*
      if( ( uexp+emin ).gt.( -lexp-emin ) ) then
         expsum = 2*lexp
      else
         expsum = 2*uexp
      end if
*
*     expsum is the exponent range, approximately equal to
*     emax - emin + 1 .
*
      emax = expsum + emin - 1
      nbits = 1 + exbits + p
*
*     nbits is the total number of bits needed to store a
*     floating-point number.
*
      if( ( mod( nbits, 2 ).eq.1 ) .and. ( beta.eq.2 ) ) then
*
*        either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. cray machines) or the mantissa has an implicit bit,
*        (e.g. ieee machines, dec vax machines), which is perhaps the
*        most likely. we have to assume the last alternative.
*        if this is true, then we need to reduce emax by one because
*        there must be some way of representing zero in an implicit-bit
*        system. on machines like cray, we are reducing emax by one
*        unnecessarily.
*
         emax = emax - 1
      end if
*
      if( ieee ) then
*
*        assume we are on an ieee machine which reserves one exponent
*        for infinity and nan.
*
         emax = emax - 1
      end if
*
*     now create rmax, the largest machine number, which should
*     be equal to (1.0 - beta**(-p)) * beta**emax .
*
*     first compute 1.0 - beta**(-p), being careful that the
*     result is less than 1.0 .
*
      recbas = one / beta
      z = beta - one
      y = zero
      do 20 i = 1, p
         z = z*recbas
         if( y.lt.one )
     $      oldy = y
         y = dlamc3( y, z )
   20 continue
      if( y.ge.one )
     $   y = oldy
*
*     now multiply by beta**emax to get rmax.
*
      do 30 i = 1, emax
         y = dlamc3( y*beta, zero )
   30 continue
*
      rmax = y
      return
*
*     end of dlamc5
*
      end


      logical          function lsame( ca, cb )
*
*  -- lapack auxiliary routine (version 3.0) --
*     univ. of tennessee, univ. of california berkeley, nag ltd.,
*     courant institute, argonne national lab, and rice university
*     september 30, 1994
*
*     .. scalar arguments ..
      character          ca, cb
*     ..
*
*  purpose
*  =======
*
*  lsame returns .true. if ca is the same letter as cb regardless of
*  case.
*
*  arguments
*  =========
*
*  ca      (input) character*1
*  cb      (input) character*1
*          ca and cb specify the single characters to be compared.
*
* =====================================================================
*
*     .. intrinsic functions ..
      intrinsic          ichar
*     ..
*     .. local scalars ..
      integer            inta, intb, zcode
*     ..
*     .. executable statements ..
*
*     test if the characters are equal
*
      lsame = ca.eq.cb
      if( lsame )
     $   return
*
*     now test for equivalence if both characters are alphabetic.
*
      zcode = ichar( 'z' )
*
*     use 'z' rather than 'a' so that ascii can be detected on prime
*     machines, on which ichar returns a value with bit 8 set.
*     ichar('a') on prime machines returns 193 which is the same as
*     ichar('a') on an ebcdic machine.
*
      inta = ichar( ca )
      intb = ichar( cb )
*
      if( zcode.eq.90 .or. zcode.eq.122 ) then
*
*        ascii is assumed - zcode is the ascii code of either lower or
*        upper case 'z'.
*
         if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
         if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
*
      else if( zcode.eq.233 .or. zcode.eq.169 ) then
*
*        ebcdic is assumed - zcode is the ebcdic code of either lower or
*        upper case 'z'.
*
         if( inta.ge.129 .and. inta.le.137 .or.
     $       inta.ge.145 .and. inta.le.153 .or.
     $       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
         if( intb.ge.129 .and. intb.le.137 .or.
     $       intb.ge.145 .and. intb.le.153 .or.
     $       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
*
      else if( zcode.eq.218 .or. zcode.eq.250 ) then
*
*        ascii is assumed, on prime machines - zcode is the ascii code
*        plus 128 of either lower or upper case 'z'.
*
         if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
         if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
      end if
      lsame = inta.eq.intb
*
*     return
*
*     end of lsame
*
      end
