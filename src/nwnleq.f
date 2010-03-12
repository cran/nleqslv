
      subroutine nwnleq(x0,n,scalex,maxit,
     *                  jacflg,xtol,ftol,btol,method,global,xscalm,
     *                  stepmx,dlt,sigma,rwork,lrwork,
     *                  rcdwrk,icdwrk,qrwork,qrwsiz,fjac,fvec,outopt,xp,
     *                  fp,gp,njcnt,nfcnt,termcd)

      integer n,jacflg,maxit,njcnt,nfcnt,termcd,method
      integer global,xscalm,lrwork,qrwsiz
      integer outopt(*)
      double precision  xtol,ftol,btol,stepmx,dlt,sigma
      double precision  xp(*),fp(*),gp(*),x0(*)
      double precision  rwork(*),rcdwrk(*),qrwork(*)
      double precision  scalex(*)         
      integer           icdwrk(*)
      external fjac,fvec

c-------------------------------------------------------------------------
c
c     Solves systems of nonlinear equations using the Newton / Broyden
c     method with a global strategy either linesearch or double dogleg
c
c     In       x0      Real(*)         starting vector for x
c     In       n       Integer         dimension of problem
c     Inout    scalex  Real(*)         scaling factors x()
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer         jacobian flag 
c                                        0 numeric
c                                        1 analytical (user supplied)
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    method  Integer         method to use 
c                                        0 Newton
c                                        1 Broyden
c     Inout    global  Integer         global strategy to use
c                                        0 quadratic linesearch
c                                        1 geometric linesearch
c                                        2 double dogleg
c                                        3 powell dogleg
c     In       xscalm  Integer         scaling method                              
c                                        0 scale fixed and supplied by user
c                                        1 for scale from jac. columns a la Minpack
c     Inout    stepmx  Real            maximum stepsize
c     Inout    dlt     Real            trust region radius 
c                                        > 0.0 or special value for initial value
c                                        -1.0  ==> use min(Cauchy length, stepmx)
c                                        -2.0  ==> use min(Newton length, stepmx)
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Out      rwork   Real(*)         real workspace (9n+2n^2)
c     In       lrwork  Integer         size real workspace
c     In       rcdwrk  Real(*)         workspace for Dtrcon (3n)
c     In       icdwrk  Integer(*)      workspace for Dtrcon (n) 
c     In       qrwork  Real(*)         workspace for Lapack QR routines (call nwqmem)
c     In       qrwsiz  Integer         size of qrwork
c     In       fjac    Name            optional name of routine to calculate
c                                      analytic jacobian
c     In       fvec    Name            name of routine to calculate f(x)
c     In       outopt  Integer(*)      output options
c                                       outopt(1)
c                                         0 no output
c                                         1 output an iteration report
c                                       outopt(2)
c                                         0 do not check any analytical jacobian
c                                         1 check analytical jacobian if supplied
c     Out      xp      Real(*)         final values for x()
c     Out      fp      Real(*)         final values for f(x)
c     Out      gp      Real(*)         gradient of f() at xp()
c     Out      njcnt   Integer         number of jacobian evaluations
c     Out      nfcnt   Integer         number of function evaluations
c     Out      termcd  Integer         termination code
c                                       > 0 process terminated                                  
c                                             1  function criterion near zero
c                                             2  no better point found
c                                             3  x-values within tolerance
c                                             4  iteration limit exceeded
c
c                                       < 0 invalid input parameters
c                                            -1  n not positive
c                                            -2  insufficient workspace rwork
c                                            -3  cannot check analytical jacobian (not supplied)
c
c    The subroutine fvec must be declared as
c
c!        subroutine fvec(x,f,n,flag)
c         double precision x(*), f(*)
c         integer  n, flag
c                   
c         x() are the x values for which to calculate the function values f(*)
c         The dimension of these vectors is n
c         The flag argument is set to
c            0  for calculation of function values
c           >0  indicating that jacobian column <flag> is being computed
c               so that fvec can abort.
c 
c    The subroutine fjac must be declared as
c
c!        subroutine mkjac(rjac,ldr,x,n)    
c         integer ldr
c         double precision rjac(ldr,*), x(*)
c         integer  n
c
c         The routine calculates the jacobian in point x(*) of the 
c         function. If any illegal values are encountered during
c         calculation of the jacobian it is the responsibility of
c         the routine to quit.

c-------------------------------------------------------------------------

      double precision epsm

c     check input parameters

      call nwpchk(n,lrwork,xtol,ftol,btol,maxit,
     *            jacflg,method,global,stepmx,dlt,sigma,
     *            epsm,outopt,scalex,termcd)
      if(termcd .lt. 0) then
         return
      endif
      
c     first argument of nwsolv/brsolv is leading dimension of rjac in those routines
c     should be at least n

      if( method .eq. 0 ) then

         call nwsolv(n,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,global,xscalm,
     *               stepmx,dlt,sigma,
     *               rwork(1+9*n),
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,
     *               epsm,fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,termcd)

      elseif( method .eq. 1 ) then

         call brsolv(n,x0,n,scalex,maxit,jacflg,
     *               xtol,ftol,btol,global,xscalm,
     *               stepmx,dlt,sigma,
     *               rwork(1+9*n),
     *               rwork(1    ),rwork(1+  n),
     *               rwork(1+2*n),rwork(1+3*n),
     *               rwork(1+4*n),rwork(1+5*n),
     *               rwork(1+6*n),rwork(1+7*n),
     *               rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,
     *               epsm,fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,termcd)

      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine nwqmem(n,wrksiz)
      integer n, wrksiz

c-------------------------------------------------------------------------
c     Query the size of the double precision work array required
c     for optimal operation of the Lapack QR routines
c-------------------------------------------------------------------------

      double precision A(1), work(1)
      integer lwork, info
      
      lwork = -1
      call dgeqrf(n,n,A,n,work,work,lwork,info)
      if( info .ne. 0 ) then
          wrksiz = -1
      else
          wrksiz = int(work(1))
      endif
      
      return
      end
      
c-----------------------------------------------------------------------

      subroutine nwpchk(n,lrwk,
     *                  xtol,ftol,btol,maxit,jacflg,method,global,
     *                  stepmx,dlt,sigma,epsm,outopt,
     *                  scalex,termcd)

      integer n,lrwk,jacflg
      integer method,global,maxit,termcd            
      integer outopt(*)
      double precision  xtol,ftol,btol,stepmx,dlt,sigma,epsm
      double precision  scalex(*)

c-------------------------------------------------------------------------
c
c     Check input arguments for consistency and modify if needed/harmless
c
c     Arguments
c
c     In       n       Integer         dimension of problem
c     In       lrwk    Integer         size real workspace
c     Inout    xtol    Real            x tolerance
c     Inout    ftol    Real            f tolerance
c     Inout    btol    Real            x tolerance for backtracking
c     Inout    maxit   Integer         maximum number iterations
c     Inout    jacflg  Integer         jacobian flag
c     Inout    method  Integer         method to use (Newton/Broyden)
c     Inout    global  Integer         global strategy to use
c     Inout    stepmx  Real            maximum stepsize
c     Inout    dlt     Real            trust region radius
c     Inout    sigma   Real            reduction factor geometric linesearch
c     Out      epsm                    machine precision
c     Inout    scalex  Real(*)         scaling factors x()
c     Out      termcd  Integer         termination code (<0 on errors)
c
c-------------------------------------------------------------------------

      integer i,len
      double precision epsmch

      double precision Rzero, Rone, Rtwo, Rthree
      parameter(Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0)

      double precision Rthous
      parameter(Rthous = 1000.0d0)

      double precision Rhalf
      parameter(Rhalf = 0.5d0)

c     check that parameters only take on acceptable values
c     if not, set them to default values

c     initialize termcd to all ok

      termcd = 0

c     compute machine precision

      epsm = epsmch()

c     check dimensions of the problem

      if(n .le. 0) then
         termcd = -1
         return
      endif

c     check dimensions of workspace arrays

      len = 9*n+2*n*n
      if(lrwk .lt. len) then
         termcd = -2
         return
      endif

c     check jacflg, method, and global

      if(jacflg .ne. 1) jacflg = 0

      if(method .lt. 0 .or. method .gt. 1) method = 0

      if(global .lt. 0 .or. global .gt. 3) global = 0

c     set outopt to correct values

      if(outopt(1) .ne. 0 ) then
         outopt(1) = 1
      endif

      if(outopt(2) .ne. 0 ) then
         outopt(2) = 1
      endif
      
c     check scale matrices

      do i = 1,n
         if(scalex(i) .lt. Rzero) scalex(i) = -scalex(i)
         if(scalex(i) .eq. Rzero) scalex(i) = Rone
      enddo

c     check step and function tolerances

      if(xtol .lt. Rzero) then
         xtol = epsm**(Rtwo/Rthree)
      endif

      if(ftol .lt. Rzero) then
         ftol = epsm**(Rtwo/Rthree)
      endif

      if( btol .lt. xtol ) btol = xtol

c     check reduction in geometric linesearch

      if( sigma .le. Rzero .or. sigma .ge. Rone ) then
         sigma = Rhalf
      endif

c     check iteration limit

      if(maxit .le. 0) then
         maxit = 150
      endif

c     check stepmx and dlt

      if(stepmx .lt. Rzero) stepmx = Rthous

      if(dlt .le. Rzero) then
         if( dlt .ne. -Rtwo ) then
            dlt = -Rone
         endif
      elseif(dlt .gt. stepmx) then
         dlt = stepmx
      endif

      return
      end
