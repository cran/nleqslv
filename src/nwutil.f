
      subroutine nwtcvg(xplus,fplus,xc,scalex,xtol,retcd,ftol,iter,
     *                  maxit,n,termcd)

      integer n,iter,maxit,termcd,retcd
      double precision  xtol,ftol
      double precision  xplus(*),fplus(*),xc(*), scalex(*)

c-------------------------------------------------------------------------
c
c     Decide whether to terminate the nonlinear algorithm
c
c     Arguments
c
c     In       xplus   Real(*)         new x values
c     In       fplus   Real(*)         new f values
c     In       xc      Real(*)         current x values
c     In       scalex  Real(*)         scaling factors x
c     In       xtol    Real            stepsize tolerance
c     In       retcd   Integer         return code from global search routines
c     In       ftol    Real            function tolerance
c     In       iter    Integer         iteration number
c     In       maxit   Integer         maximum number of iterations allowed
c     In       n       Integer         size of x and f
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
c
c-------------------------------------------------------------------------

      double precision  fmax,rsx, nuxnrm
      integer idamax

c     check whether function values are within tolerance

      fmax = abs(fplus(idamax(n,fplus,1)))
      if( fmax .le. ftol) then
         termcd = 1
         goto 1000
      endif

      if(iter .eq. 0) goto 1000

      if(retcd .eq. 1) then
         termcd = 3
         goto 1000
      endif

c     check whether relative step length is within tolerance
c     Dennis Schnabel Algorithm A7.2.3

      rsx = nuxnrm(n, xplus, xc, scalex)
      if(rsx .le. xtol) then
        termcd = 2
        goto 1000
      endif

c     check iteration limit

      if(iter .ge. maxit) then
         termcd = 4
      endif

 1000 continue
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
      double precision  ndigit,rnoise,stepsz,xtmpj,dinf
      double precision  tol
      integer idamax

      integer MAXERR
      parameter(MAXERR=10)

      double precision Rquart, Rone, Rten
      parameter(Rquart=0.25d0, Rone=1.0d0, Rten=10.0d0)

      termcd = 0

c     compute the finite difference jacobian and check it against
c     the analytic one

      ndigit = -log10(epsm)
      rnoise = max(Rten**(-ndigit),epsm)
      rnoise = sqrt(rnoise)
      tol    = epsm**Rquart

      errcnt = 0

      do 40 j = 1,n
         stepsz = rnoise*max(abs(xc(j)),Rone)
         xtmpj = xc(j)
         xc(j) = xtmpj+stepsz
         call fvec(xc,fz,n,0)
         xc(j) = xtmpj

         do 10 i = 1,n
            wa(i) = (fz(i)-fc(i))/stepsz
 10      continue

         dinf = abs(wa(idamax(n,wa,1)))

         do 30 i = 1,n
            if(abs(A(i,j)-wa(i)).gt.tol*dinf) then
               errcnt = errcnt + 1
               if( errcnt .gt. MAXERR ) then
                     goto 50
               endif
               call nwckot(i,j,A(i,j),wa(i))
            endif
 30      continue

 40   continue

 50   continue

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

      do 20 j = 1,n
         h  = p + p * abs(xc(j))

c        or as alternative h  = p * max(Rone, abs(xc(j)))

         xcj   = xc(j)
         xc(j) = xcj + h

c        avoid (small) rounding errors
c        h = xcj - xc(j) but not here to avoid clever optimizers

         h = rnudif(xc(j), xcj)
         call fvec(xc,fz,n,j)
         xc(j) = xcj
         do 10 i = 1,n
            rjac(i,j) = (fz(i)-fc(i)) / h
 10      continue
 20   continue

      return
      end

c-----------------------------------------------------------------------

      double precision function nudnrm(n, d, x, scalex)
      integer n
      double precision  d(*), x(*), scalex(*)

c-------------------------------------------------------------------------
c
c     calculate  max( abs(d[*]) / max(x[*],1/scalex[*]) )
c
c     Arguments
c
c     In   n        Integer       number of elements in d() and x()
c     In   d        Real(*)       vector d
c     In   x        Real(*)       vector x
c     In   scalex   Real(*)       scaling vector sx
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do 10 i=1,n
         t = max( t, abs(d(i)) / max(abs(x(i)),Rone/scalex(i)) )
  10  continue

      nudnrm = t

      return
      end

c-----------------------------------------------------------------------

      double precision function nuxnrm(n, xn, xc, scalex)
      integer n
      double precision  xn(*), xc(*), scalex(*)

c-------------------------------------------------------------------------
c
c     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1/scalex[*]) )
c
c     Arguments
c
c     In   n        Integer       number of elements in xn() and xc()
c     In   xn       Real(*)       vector xn
c     In   xc       Real(*)       vector xc
c     In   scalex   Real(*)       scaling factors x(*)
c
c-------------------------------------------------------------------------

      integer i
      double precision  t

      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

      t = Rzero
      do 10 i=1,n
         t = max( t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone/scalex(i)) )
  10  continue

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

      subroutine compmu(r,ldr,n,epsm,mu,y)

      integer ldr,n
      double precision r(ldr,*),epsm,mu,y(*)

c-------------------------------------------------------------------------
c
c     Compute a small perturbation mu for the (almost) singular matrix R.
c     mu is used in the computation of the levenberg-marquardt step.
c
c     Arguments
c
c     In       R       Real(ldr,*)     upper triangular matrix from QR
c     In       ldr     Integer         leading dimension of R
c     In       n       Integer         column dimension of R
c     In       epsm    Real            machine precision
c     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
c                                      * n * epsm * 100) designed to make
c                                        trans(R)*R + mu * I not singular
c     Wk       y       Real(*)         workspace for dlange
c
c-------------------------------------------------------------------------

      double precision  aifnrm,al1nrm
      double precision dlantr

      double precision Rhund
      parameter(Rhund=100d0)

c     get the infinity norm of R
c     get the l1 norm of R

      aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
      al1nrm = dlantr('1','U','N',n,n,r,ldr,y)

c     compute mu

      mu = sqrt(aifnrm*al1nrm*n*epsm*Rhund)

      return
      end

c-----------------------------------------------------------------------

      subroutine cndjac(n,r,ldr,epsm,rcond,y,rcdwrk,icdwrk,ierr,mu)
      integer n,ldr,icdwrk(*),ierr
      double precision epsm,rcond,mu,r(ldr,*),y(*),rcdwrk(*)

c---------------------------------------------------------------------
c
c     Check r for singularity and/or ill conditioning
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       r       Real(ldr,*)     upper triangular R from QR decomposition
c     In       ldr     Integer         leading dimension of rjac
c     In       epsm    Real            machine precision
c     Out      rcond   Real            inverse condition  of r
c     Wk       y       Real(*)         workspace
c     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
c     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
c     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
c                                      1 indicating Jacobian ill-conditioned
c                                      2 indicating Jacobian completely singular
c     Out      mu      Real            0 if ierr == 0
c                                      small positive number when ierr > 0
c                                      to make trans(R)*R+mu*I non singular
c
c---------------------------------------------------------------------

      integer i,info
      logical rsing
      double precision Rzero,R2d3
      parameter(Rzero=0.0d0, R2d3=2.0d0/3.0d0)

      mu = Rzero
      ierr = 0

      rsing = .false.
      do 10 i=1,n
         if( r(i,i) .eq. Rzero ) then
             rsing = .true.
         endif
   10 continue

      if( rsing ) then
         ierr = 2
         rcond = Rzero
      else
         call dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         if( rcond .lt. epsm**R2d3 ) then
             ierr = 1
         endif
      endif

      if( ierr .gt. 0 ) then
         call compmu(r,ldr,n,epsm,mu,y)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfjac(x,f,fq,n,epsm,jacflg,fvec,mkjac,rjac,ldr)

      integer ldr,n,jacflg
      double precision  epsm
      double precision  x(*),f(*)
      double precision  rjac(ldr,*),fq(*)
      external fvec,mkjac

c-------------------------------------------------------------------------
c
c     Calculate and scales the jacobian  matrix
c
c     Arguments
c
c     In       x       Real(*)         (scaled) current x values
c     In       f       Real(*)         function values f(x)
c     Wk       fq      Real(*)         (internal) workspace
c     In       n       Integer         size of x and f
c     In       epsm    Real            machine precision
c     In       jacflg  Integer         indicates how to compute jacobian
c                                       0  numeric
c                                       1  analytic
c     In       fvec    Name            name of routine to evaluate f()
c     In       mkjac   Name            name of routine to evaluate jacobian
c     Out      rjac    Real(ldr,*)     jacobian matrix
c     In       ldr     Integer         leading dimension of rjac
c
c-------------------------------------------------------------------------

      integer i,j
      double precision  t
      double precision  dnrm2

c     compute the finite difference or analytic jacobian at x

      if(jacflg .eq. 0) then
         call fdjac(x,f,n,epsm,fvec,fq,rjac,ldr)
      else
         call mkjac(rjac,ldr,x,n)
      endif

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
c     Out      stepsiz Real       stepsize   
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

      do 10 i = 1,n
         x(i) = sx(i) * x(i)
 10   continue

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

      do 10 i = 1,n
         x(i) = x(i) / sx(i)
 10   continue

      return
      end

c-----------------------------------------------------------------------

      subroutine nwfvec(x,n,fvec,f,fnorm)

      integer n
      double precision  x(*),f(*),fnorm
      external fvec

c-------------------------------------------------------------------------
c
c     Evaluate the function at current iterate x and scale its value
c
c     Arguments
c
c     In       x       Real(*)         x
c     In       n       Integer         size of x
c     In       fvec    Name            name of routine to calculate f(x)
c     Out      f       Real(*)         f(x)
c     Out      fnorm   Real            .5*||f(x)||**2
c
c-------------------------------------------------------------------------

      double precision dnrm2

      double precision Rhalf
      parameter(Rhalf=0.5d0)

      call fvec(x,f,n,0)

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
