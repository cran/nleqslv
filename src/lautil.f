
c-----------------------------------------------------------------------------

      subroutine liqrfa(a, lda, n, tau, work, wsiz, info)
      integer  lda, n, wsiz, info
      double precision  a(lda,*), tau(*), work(*)

c-------------------------------------------------------------
c
c     QR decomposition of A(n,n)
c
c     Arguments
c
c      Inout A        Real(Lda,n)    Matrix to transform.
c      In    lda      Integer        Leading dimension of A
c      In    n        Integer        number of rows/cols A
c      Out   tau      Real(n)        Information for recovering
c      Out   work     Real(*)        workspace
c      In    wsiz     Integer        size of work()

c     Lapack blocked QR (much faster for larger n)
c
c-------------------------------------------------------------

      call dgeqrf(n,n,A,lda,tau,work,wsiz,info)

      return
      end

c=============================================================

      subroutine liqsiz(n,wrksiz)
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

c=============================================================

      subroutine liqrqt(a, lda, n, tau, qty, work, wsiz, info)
      integer lda, n, wsiz, info
      double precision a(lda,*), tau(*), qty(*), work(*)

c-------------------------------------------------------------
c      Arguments
c
c      In    A     Real(Lda, n)    QR decomposition
c      In    Lda   Integer         Leading dimension A
c      In    n     Integer         Number of rows/columns in A
c      In    tau   Integer         Householder constants from QR
c      Inout qty   Real(n)         On input, vector y
c                                  On output, trans(Q)*y
c      Out   work  Real(*)         workspace
c      In    wsiz  Integer         size of work()
c
c     Liqrqt calculates trans(Q)*y from the QR decomposition
c
c     Lapack blocked
c-------------------------------------------------------------

      call dormqr('L','T',n,1,n,A,lda,tau,qty,n,work,wsiz,info)

      return
      end

c=============================================================

      subroutine liqrqq(q,ldq,tau,n,work,wsiz,info)
      integer n, ldq, wsiz, info
      double precision  q(ldq,*),tau(*),work(*)

c     Arguments
c
c     Inout  Q     Real(ldq,*)     On input, QR decomposition
c                                    On output, the full Q
c     In     ldq   Integer         leading dimension of Q
c     In     tau   Real(n)         Householder constants of
c                                     the QR decomposition
c     In     n     Integer         number of rows/columns in Q
c     Out    work  Real(n)         workspace of length n
c     In     wsiz  Integer         size of work()
c
c     Generate Q from QR decomposition Liqrfa (dgeqr2)
c
c     Lapack blocked
c-------------------------------------------------------------

      call dorgqr(n,n,n,q,ldq,tau,work,wsiz,info)

      return
      end

c-----------------------------------------------------------------------------

      subroutine nuzero(n,x)
      integer n
      double precision x(*)

c     Parameters:
c
c     In    n        Integer           Number of elements.
c     In    x        Real(*)           Vector of reals.
c
c     Description:
c
c     Nuzero sets all elements of x to 0.
c     Does nothing when n <= 0

      double precision Rzero
      parameter(Rzero=0.0d0)

      integer i

      do i=1,n
         x(i) = Rzero
      enddo

      return
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
c     Description
c
c     Nuvgiv calculates the givens rotator
c
c             |  c   s |
c         G = |        |
c             | -s   c |
c
c     with  c*c+s*s=1
c
c     for which G * | x | = | t |
c                   | y |   | 0 |
c
c     resulting in
c
c            c * x + s * y = t
c           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
c
c     Use Lapack dlartg routine
c     return c and s and the modified x and y
c     This differs from dlartg which does not modify input arguments.
c     See http://www.netlib.org/lapack/explore-html/dd/d24/dlartg_8f_source.html
c     c * x + s * y may differ from t with machine precision

      double precision t

      double precision Rzero
      parameter(Rzero=0.0d0)

      call dlartg(x,y,c,s,t)
      x = t
      y = Rzero
      return
      end
