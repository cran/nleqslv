c ----------------------------------------------------------------------

      subroutine nuvset(n,val,x,incx)
      integer n,incx
      double precision x(*),val

c  Parameters:
c
c      In    n        Integer           Number of elements.
c      In    val      Real              constant
c      In    x        Real(*)           Vector of reals.
c      In    incx     Integer           Steplength in x.
c
c  Description:
c
c      Subroutine Nuvset sets all elements of x to val.
c
c      Does nothing when n <= 0 or Incx <= 0

      integer ns,m,i

      if(n .le. 0 .or. incx .le. 0) goto 100

      if( incx .eq. 1 ) then

c     **** code for increment equal to 1
c     **** clean-up loop so remaining vector length is a multiple of 7.

         m = mod(n,7)
         do 10 i = 1,m
            x(i) = val
   10    continue
         do 20 i = m+1,n,7
            x(i)     = val
            x(i + 1) = val
            x(i + 2) = val
            x(i + 3) = val
            x(i + 4) = val
            x(i + 5) = val
            x(i + 6) = val
   20    continue

      else

c     **** code for positive, nonunit increment
         ns = n * incx
         do 30 i=1,ns,incx
            x(i) = val
   30    continue

      endif

  100 return
      end

c ----------------------------------------------------------------------

      subroutine nuvgiv(x,y,c,s)
      double precision x,y,c,s

c     Parameters
c
c     Inout   x     Real       x input / c*x+s*y on output
c     Inout   y     Real       y input / 0       on output
c     Out     c     Real       c of tranformation (cosine)
c     Out     s     Real       s of tranformation (  sine)
c
c
c     Description
c
c     Nuvgiv calculates the givens rotator
c
c
c             |  c   s |
c         G = |        |
c             | -s   c |
c
c     with  c*c+s*s=1
c
c     for which G * | x | = | z |
c                   | y |   | 0 |
c
c     then we have
c
c            c * x + s * y = z
c           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
c            
c     Use Lapack dlartg routine  
c     and return c and s modified x and y

      double precision t

      double precision Rzero
      parameter(Rzero=0.0d0)
      
      call dlartg(x,y,c,s,t)  
      x = t
      y = Rzero
      return
      end
