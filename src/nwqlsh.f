
      subroutine nwqlsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,fvec,
     *                  xp,fp,fpnorm,mxtake,retcd,gcnt,priter,iter)

      integer n,retcd,gcnt
      double precision  stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(*),fp(*)
      double precision  scalex(*)
      logical mxtake
      external fvec

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using a safeguarded quadratic line search
c     along the newton direction
c
c     Arguments
c
c     In       n       Integer          dimension of problem
c     In       xc      Real(*)          current iterate
c     In       fcnorm  Real             0.5 * || f(xc) ||**2
c     In       d       Real(*)          newton direction
c     In       g       Real(*)          gradient at current iterate
c     In       stepmx  Real             maximum stepsize
c     In       xtol    Real             relative step size at which
c                                       successive iterates are considered
c                                       close enough to terminate algorithm
c     In       scalex  Real(*)          diagonal scaling matrix for x()
c     In       fvec    Name             name of routine to calculate f()
c     In       xp      Real(*)          new x()
c     In       fp      Real(*)          new f(x)
c     In       fpnorm  Real             .5*||fp||**2
c
c     Out      mxtake  Logical          .true. if maximum step taken
c                                       else .false.
c
c     Out      retcd   Integer          return code
c                                         0 new satisfactory x() found
c                                         1 no  satisfactory x() found
c                                           sufficiently distinct from xc()
c
c     Out      gcnt    Integer          number of steps taken
c     In       priter  Integer           >0 unit if intermediate steps to be printed
c                                        -1 if no printing
c
c-------------------------------------------------------------------------

      integer i
      double precision  alpha,slope,rsclen,oarg(4)
      double precision  lambda,lamhi,lamlo,t
      double precision  ddot,dnrm2, nudnrm
      double precision  dlen
      logical scstep
      
      integer idamax
      
      parameter (alpha = 1.0d-4)

      double precision Rone, Rtwo, Rten
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0)

c     safeguard initial step size
c     use xp temporarily

      call dcopy(n,d,1,xp,1)
      call vscal(n,xp,scalex)
      dlen = dnrm2(n,xp,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
          scstep = .true.
      else
          lamhi  = Rone
          scstep = .false.
      endif

c     compute slope  =  g-trans * d
c     if the jacobian or an approximation is not singular or
c     ill conditioned then slope = -2*fcnorm
c     but otherwise the parameter mu used for the quasi pseudo inverse
c     destroys this equality. So calculate slope from its definition

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = nudnrm(n,d,xc,scalex)
      lamlo  = xtol / rsclen

c     initialization of retcd, mxtake and lambda (linesearch length)

      retcd  = 0
      mxtake = .false.
      lambda = lamhi
      gcnt   = 0

  50  continue

c     compute the next iterate xp

      do 60 i = 1,n
         xp(i) = xc(i) + lambda*d(i)
  60  continue

c     evaluate functions and the objective function at xp

      call nwfvec(xp,n,fvec,fp,fpnorm)
      gcnt = gcnt + 1

      if( priter .gt. 0) then
         oarg(1) = lambda
         oarg(2) = fcnorm + alpha * lambda * slope
         oarg(3) = fpnorm
         oarg(4) = abs(fp(idamax(n,fp,1)))
         call nwlsot(iter,1,oarg)
      endif

c     test whether the standard step produces enough decrease
c     the objective function.
c     If not update lambda and compute a new next iterate

      if( fpnorm .gt. (fcnorm + alpha * lambda * slope) ) then
         t = ((-lambda**2)*slope/Rtwo)/(fpnorm-fcnorm-lambda*slope)
         lambda  = max(lambda / Rten , t)
         if(lambda .lt. lamlo) then
            retcd = 1
            goto 100
         endif
         goto 50
      endif

      if( lambda .eq. lamhi .and. scstep ) then
         mxtake = .true.
      endif

 100  continue

      return
      end
