
      subroutine nwclsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,fvec,
     *                  xp,fp,fpnorm,xw,mxtake,retcd,gcnt,priter,iter)

      integer n,retcd,gcnt
      double precision  stepmx,xtol,fcnorm,fpnorm
      double precision  xc(*)
      double precision  d(*),g(*),xp(*),fp(*),xw(*)
      double precision  scalex(*)
      logical mxtake
      external fvec

      integer priter,iter

c-------------------------------------------------------------------------
c
c     Find a next acceptable iterate using a safeguarded cubic line search
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
c     Out      xw      Real(*)           workspace for unscaling x(*)
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
      double precision  ddot,dnrm2, nudnrm, ftarg
      double precision  dlen
      logical scstep
      double precision a, b, disc, fpt, fpt0, fpnorm0, lambda0
      integer idamax
      logical firstback
      
      parameter (alpha = 1.0d-4)

      double precision Rone, Rtwo, Rthree, Rten, Rzero
      parameter(Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0, Rzero=0.0D0)
      parameter(Rthree=3.0d0)
      
      double precision dlamch
      
c     safeguard initial step size

      dlen = dnrm2(n,d,1)
      if( dlen .gt. stepmx ) then
          lamhi  = stepmx / dlen
          scstep = .true.
      else
          lamhi  = Rone
          scstep = .false.
      endif

c     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

c     compute the smallest value allowable for the damping
c     parameter lambda ==> lamlo

      rsclen = nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

c     initialization of retcd, mxtake and lambda (linesearch length)

      retcd  = 2
      mxtake = .false.
      lambda = lamhi
      gcnt   = 0
      firstback = .true.
      
      do while( retcd .eq. 2 )

c        compute next x

         do i=1,n
            xp(i) = xc(i) + lambda*d(i)
         enddo

c        evaluate functions and the objective function at xp

         call nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

         if( priter .gt. 0) then
            oarg(1) = lambda
            oarg(2) = ftarg
            oarg(3) = fpnorm
            oarg(4) = abs(fp(idamax(n,fp,1)))
            call nwlsot(iter,1,oarg)
         endif

c        first is quadratic
c        test whether the standard step produces enough decrease
c        of the objective function.
c        If not update lambda and compute a new next iterate

         if( fpnorm .le. ftarg ) then
             retcd = 0
         else
            if( fpnorm .gt. lamlo**2 * sqrt(dlamch('O')) ) then
c               safety against overflow in what folows (use lamlo**2 for safety)
                lambda = lambda/Rten
                firstback = .true.
            else
                if( firstback ) then
                   t = ((-lambda**2)*slope/Rtwo)/
     *                        (fpnorm-fcnorm-lambda*slope)
                   firstback = .false.
                else
                   fpt  = fpnorm -fcnorm - lambda*slope
                   fpt0 = fpnorm0-fcnorm - lambda0*slope
                   a = fpt/lambda**2 - fpt0/lambda0**2
                   b = -lambda0*fpt/lambda**2 + lambda*fpt0/lambda0**2
                   a = a /(lambda - lambda0)
                   b = b /(lambda - lambda0)
                   disc = b**2 - Rthree * a * slope  
                   if( abs(a) .le. dlamch('E') ) then
                       t = -slope/(2*b)
                   elseif(disc .gt. b**2) then
c                      a single positive solution                       
                       t = (-b + sign(Rone,a)*sqrt(disc))/(Rthree*a)
                   else
c                      both roots > 0, left one is minimum
                       t = (-b - sign(Rone,a)*sqrt(disc))/(Rthree*a)
                   endif
                   t = min(t, .5*lambda)
                endif
                lambda0 = lambda
                fpnorm0 = fpnorm
                lambda = max(t,lambda/Rten)
                if(lambda .lt. lamlo) then
                   retcd = 1
                endif
            endif
         endif 
      enddo

      if( lambda .eq. lamhi .and. scstep ) then
         mxtake = .true.
      endif

      return
      end
