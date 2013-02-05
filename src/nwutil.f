
      subroutine nwtcvg(xplus,fplus,xc,xtol,retcd,ftol,iter,
     *                  maxit,n,ierr,termcd)

      integer n,iter,maxit,ierr,termcd,retcd
      double precision  xtol,ftol
      double precision  xplus(*),fplus(*),xc(*)

c-------------------------------------------------------------------------
c
c     Decide whether to terminate the nonlinear algorithm
c
c     Arguments
c
c     In       xplus   Real(*)         new x values
c     In       fplus   Real(*)         new f values
c     In       xc      Real(*)         current x values
c     In       xtol    Real            stepsize tolerance
c     In       retcd   Integer         return code from global search routines
c     In       ftol    Real            function tolerance
c     In       iter    Integer         iteration number
c     In       maxit   Integer         maximum number of iterations allowed
c     In       n       Integer         size of x and f
c     In       ierr    Integer         return code of cndjac (condition estimation)
c
c     Out      termcd  Integer         termination code
c                                        0 no termination criterion satisfied
c                                          ==> continue iterating
c                                        1 norm of scaled function value is
c                                          less than ftol
c                                        2 scaled distance between last
c                                          two steps less than xtol
c                                        3 unsuccessful global strategy
c                                          ==> cannot find a better point
c                                        4 iteration limit exceeded
c                                        5 Jacobian too ill-conditioned
c                                        6 Jacobian singular
c
c-------------------------------------------------------------------------

      double precision  fmax,rsx, nuxnrm
      integer idamax

c     check whether function values are within tolerance

      termcd = 0

      if( ierr .ne. 0 ) then
         termcd = 4 + ierr
         return
      endif

      fmax = abs(fplus(idamax(n,fplus,1)))
      if( fmax .le. ftol) then
         termcd = 1
         return
      endif

c     initial check at start so there is no xplus
c     so only a check of function values is useful
      if(iter .eq. 0) return

      if(retcd .eq. 1) then
         termcd = 3
         return
      endif

c     check whether relative step length is within tolerance
c     Dennis Schnabel Algorithm A7.2.3

      rsx = nuxnrm(n, xplus, xc)
      if(rsx .le. xtol) then
        termcd = 2
        return
      endif

c     check iteration limit

      if(iter .ge. maxit) then
         termcd = 4
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
      double precision xc(*),fc(*),fcnorm,xp(*),fp(*),fpnorm
      integer n, gcnt, priter, iter

c-------------------------------------------------------------------------
c
c     calling routine got an error in decomposition/update of Jacobian/Broyden
c     jacobian an singular or too ill-conditioned
c     prepare return arguments
c
c     Arguments
c
c     In       n       Integer         size of x
c     In       xc      Real(*)         current (starting) x values
c     In       fc      Real(*)         function values f(xc)
c     In       fcnorm  Real            norm fc
c     Out      xp      Real(*)         final x values
c     Out      fp      Real(*)         function values f(xp)
c     Out      fpnorm  Real            final norm fp
c     Out      gcnt    Integer         # of backtracking steps (here set to 0)
c     In       priter  Integer         flag for type of output
c     In       iter    Integer         iteration counter
c
c-------------------------------------------------------------------------

      call dcopy(n,xc,1,xp,1)
      call dcopy(n,fc,1,fp,1)
      fpnorm = fcnorm
      gcnt   = 0
      if( priter .gt. 0 ) then
         call nwjerr(iter)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine chkjac(A,lda,xc,fc,n,epsm,scalex,fz,wa,fvec,termcd)

      integer lda,n,termcd
      double precision  A(lda,*),xc(*),fc(*)
      double precision  epsm,scalex(*)
      double precision  fz(*),wa(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Check the analytic jacobian against its finite difference approximation
c
c     Arguments
c
c     In       A       Real(lda,*)     analytic jacobian
c     In       lda     Integer         leading dimension of ajanal
c     In       xc      Real(*)         vector of x values
c     In       fc      Real(*)         function values f(xc)
c     In       n       Integer         size of x
c     In       epsm    Real            machine precision
c     In       scalex  Real(*)         scaling vector for x()
c     Wk       fz      Real(*)         workspace
c     Wk       wa      Real(*)         workspace
c     In       fvec    Name            name of routine to evaluate f(x)
c     Out      termcd  Integer         return code
c                                        0  analytic jacobian ok
c                                      -10  analytic jacobian NOT ok
c
c-------------------------------------------------------------------------

      integer i,j,errcnt
      double precision  ndigit,p,h,xcj,dinf
      double precision  tol
      double precision  rnudif
      integer idamax

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rone, Rten
      parameter(Rquart=0.25d0, Rone=1.0d0, Rten=10.0d0)

      termcd = 0

c     compute the finite difference jacobian and check it against
c     the analytic one

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))
      tol    = epsm**Rquart

      errcnt = 0
      call vunsc(n,xc,scalex)

      do j=1,n
         h = p + p * abs(xc(j))
         xcj   = xc(j)
         xc(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj but not here to avoid clever optimizers

         h = rnudif(xc(j), xcj)

         call fvec(xc,fz,n,j)
         xc(j) = xcj

         do i=1,n
            wa(i) = (fz(i)-fc(i))/h
         enddo

         dinf = abs(wa(idamax(n,wa,1)))

         do i=1,n
            if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
               errcnt = errcnt + 1
               if( errcnt .gt. MAXERR ) then
                  termcd = -10
                  return
               endif
               call nwckot(i,j,A(i,j),wa(i))
            endif
         enddo
      enddo

      call vscal(n,xc,scalex)

      if( errcnt .gt. 0 ) then
         termcd = -10
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine fdjac(xc,fc,n,epsm,fvec,fz,rjac,ldr)

      integer ldr,n
      double precision  epsm
      double precision  rjac(ldr,*),fz(*),xc(*),fc(*)
      external fvec

c-------------------------------------------------------------------------
c
c     Compute the finite difference jacobian at the current point xc
c
c     Arguments
c
c     In       xc      Real(*)         current point
c     In       fc      Real(*)         function values at current point
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       fvec    Name            name of routine to evaluate f(x)
c     Wk       fz      Real(*)         workspace
c     Out      rjac    Real(ldr,*)     jacobian matrix at x
c                                        entry [i,j] is derivative of
c                                        f(i) wrt to x(j)
c     In       ldr     Integer         leading dimension of rjac
c
c-------------------------------------------------------------------------

      integer i,j
      double precision  ndigit,p,h,xcj
      double precision  rnudif

      double precision Rten
      parameter(Rten=10d0)

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))

      do j=1,n
         h = p + p * abs(xc(j))

c        or as alternative h  = p * max(Rone, abs(xc(j)))

         xcj   = xc(j)
         xc(j) = xcj + h

c        avoid (small) rounding errors
c        h = xc(j) - xcj  but not here to avoid clever optimizers

         h = rnudif(xc(j), xcj)
         call fvec(xc,fz,n,j)
         xc(j) = xcj
         do i=1,n
            rjac(i,j) = (fz(i)-fc(i)) / h
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      double precision function nudnrm(n, d, x)
      integer n
      double precision  d(*), x(*)

c-------------------------------------------------------------------------
c
c     calculate  max( abs(d[*]) / max(x[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in d() and x()
c     In   d        Real(*)       vector d
c     In   x        Real(*)       vector x
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
      enddo
      nudnrm = t

      return
      end

c-----------------------------------------------------------------------

      double precision function nuxnrm(n, xn, xc)
      integer n
      double precision  xn(*), xc(*)

c-------------------------------------------------------------------------
c
c     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1) )
c
c     Arguments
c
c     In   n        Integer       number of elements in xn() and xc()
c     In   xn       Real(*)       vector xn
c     In   xc       Real(*)       vector xc
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do i=1,n
         t = max(t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone))
      enddo
      nuxnrm = t

      return
      end

c-----------------------------------------------------------------------

      double precision function rnudif(x, y)
      double precision x, y

c-------------------------------------------------------------------------
c
c     Return difference of x and y (x - y)
c
c     Arguments
c
c     In   x  Real      argument 1
c     In   y  Real      argument 2
c
c-------------------------------------------------------------------------

      rnudif = x - y
      return
      end

c-----------------------------------------------------------------------

      subroutine cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)
      integer n,ldr,icdwrk(*),ierr
      double precision cndtol,rcond,r(ldr,*),rcdwrk(*)

c---------------------------------------------------------------------
c
c     Check r for singularity and/or ill conditioning
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       r       Real(ldr,*)     upper triangular R from QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       cndtol  Real            tolerance of test for ill conditioning
c                                       when rcond <= cndtol then ierr is set to 1
c                                       cndtol should be >= machine precision
c     Out      rcond   Real            inverse condition  of r
c     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
c     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian too ill-conditioned
c                                      2 indicating Jacobian completely singular
c
c---------------------------------------------------------------------

      integer i,info
      logical rsing
      double precision Rzero,R2d3
      parameter(Rzero=0.0d0, R2d3=2.0d0/3.0d0)

      ierr = 0

      rsing = .false.
      do i=1,n
         if( r(i,i) .eq. Rzero ) then
             rsing = .true.
         endif
      enddo

      if( rsing ) then
         ierr = 2
         rcond = Rzero
      else
         call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         if( rcond .eq. Rzero ) then
             ierr = 2
         elseif( rcond .le. cndtol ) then
             ierr = 1
         endif
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfjac(x,scalex,f,fq,n,epsm,jacflg,fvec,mkjac,rjac,
     *                  ldr,xw)

      integer ldr,n,jacflg
      double precision  epsm
      double precision  x(*),f(*),scalex(*),xw(*)
      double precision  rjac(ldr,*),fq(*)
      external fvec,mkjac

c-------------------------------------------------------------------------
c
c     Calculate the jacobian  matrix
c
c     Arguments
c
c     In       x       Real(*)         (scaled) current x values
c     In       scalex  Real(*)         scaling factors x
c     In       f       Real(*)         function values f(x)
c     Wk       fq      Real(*)         (internal) workspace
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       jacflg  Integer         indicates how to compute jacobian
c                                       0  numeric
c                                       1  analytic
c     In       fvec    Name            name of routine to evaluate f()
c     In       mkjac   Name            name of routine to evaluate jacobian
c     Out      rjac    Real(ldr,*)     jacobian matrix (unscaled)
c     In       ldr     Integer         leading dimension of rjac
c     Internal xw      Real(*)         used for storing unscaled x
c
c-------------------------------------------------------------------------

c     compute the finite difference or analytic jacobian at x

      call dcopy(n,x,1,xw,1)
      call vunsc(n,xw,scalex)
      if(jacflg .eq. 0) then
         call fdjac(xw,f,n,epsm,fvec,fq,rjac,ldr)
      else
         call mkjac(rjac,ldr,xw,n)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Scale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t, Rone
      parameter(Rone=1.0d0)

      do j = 1,n
         t = Rone/scalex(j)
         call dscal(n,t,rjac(1,j),1)
      enddo 

      return
      end

c-----------------------------------------------------------------------

      subroutine nwunscjac(n,rjac,ldr,scalex)
      integer n, ldr
      double precision rjac(ldr,*), scalex(*)

c-------------------------------------------------------------------------
c
c     Unscale jacobian
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     Inout    rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     In       scalex  Real(*)         scaling factors for x
c
c-------------------------------------------------------------------------

      integer j
      double precision t

      do j = 1,n
         t = scalex(j)  
         call dscal(n,t,rjac(1,j),1)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwcpsx(n,rjac,ldr,scalex,epsm, mode)

      integer ldr,n,mode
      double precision  epsm
      double precision  scalex(*)
      double precision  rjac(ldr,*)

c-------------------------------------------------------------------------
c
c     Calculate scaling factors from the jacobian  matrix
c
c     Arguments
c
c     In       n       Integer         size of x and f
c     In       rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c     Out      scalex  Real(*)         scaling factors for x
c     In       epsm    Real            machine precision
c     In       mode    Integer         1: initialise, >1: adjust
c-------------------------------------------------------------------------

      integer k
      double precision  dnrm2

      if( mode .eq. 1 ) then
         do k=1,n
            scalex(k) = dnrm2(n,rjac(1,k),1)
            if( scalex(k) .le. epsm ) scalex(k) = 1
         enddo
      else if( mode .gt. 1 ) then
         do k=1,n
            scalex(k) = max(scalex(k),dnrm2(n,rjac(1,k),1))
         enddo
      endif
      return
      end

c-----------------------------------------------------------------------

      subroutine nwcpmt(n, x, scalex, factor, wrk, stepsiz)
      double precision x(*), scalex(*), wrk(*)
      double precision factor, stepsiz
      integer n

c-------------------------------------------------------------------------
c
c     Calculate maximum stepsize
c
c     Arguments
c
c     In       n       Integer     size of x
c     In       x       Real(*)     x-values
c     In       scalex  Real(*)     scaling factors for x
c     In       factor  Real        multiplier
c     Inout    wrk     Real(*)     workspace
c     Out      stepsiz Real        stepsize
c
c     Currently not used
c     Minpack uses this to calculate initial trust region size
c     Not (yet) used in this code because it doesn't seem to help
c     Manually setting an initial trust region size works better
c
c-------------------------------------------------------------------------

      double precision Rzero
      parameter(Rzero=0.0d0)

      double precision  dnrm2

      call dcopy(n,x,1,wrk,1)
      call vscal(n,wrk,scalex)
      stepsiz = factor * dnrm2(n,wrk,1)
      if( stepsiz .eq. Rzero ) stepsiz = factor
      return
      end

c-----------------------------------------------------------------------

      subroutine vscal(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Scale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to scale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = sx(i) * x(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine vunsc(n,x,sx)

      integer n
      double precision  x(*),sx(*)

c-------------------------------------------------------------------------
c
c     Unscale a vector x
c
c     Arguments
c
c     In       n       Integer         size of x
c     Inout    x       Real(*)         vector to unscale
c     In       sx      Real(*)         scaling vector
c
c-------------------------------------------------------------------------

      integer i

      do i = 1,n
         x(i) = x(i) / sx(i)
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfvec(x,n,scalex,fvec,f,fnorm,xw)

      integer n
      double precision  x(*),xw(*),scalex(*),f(*),fnorm
      external fvec

c-------------------------------------------------------------------------
c
c     Evaluate the function at current iterate x and scale its value
c
c     Arguments
c
c     In       x       Real(*)         x
c     In       n       Integer         size of x
c     In       scalex  Real(*)         scaling vector for x
c     In       fvec    Name            name of routine to calculate f(x)
c     Out      f       Real(*)         f(x)
c     Out      fnorm   Real            .5*||f(x)||**2
c     Internal xw      Real(*)         used for storing unscaled xc
c
c-------------------------------------------------------------------------

      double precision dnrm2

      double precision Rhalf
      parameter(Rhalf=0.5d0)

      call dcopy(n,x,1,xw,1)
      call vunsc(n,xw,scalex)
      call fvec(xw,f,n,0)

      fnorm = Rhalf * dnrm2(n,f,1)**2

      return
      end

c-----------------------------------------------------------------------

      function epsmch()

c     Return machine precision
c     Use Lapack routine

      double precision epsmch
      double precision dlamch
      external dlamch

c     dlamch('e') returns negeps (1-eps)
c     dlamch('p') returns 1+eps

      epsmch = dlamch('p')

      return
      end

c-----------------------------------------------------------------------

      function dblhuge()

c     Return largest double precision number
c     Use Lapack routine

      double precision dblhuge
      double precision dlamch
      external dlamch

c     dlamch('o') returns max double precision

      dblhuge = dlamch('o')

      return
      end
