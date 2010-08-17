      subroutine liqrup(q,ldq,n,r,ldr,u,v,wk)
      integer ldq,n,ldr
      double precision q(ldq,*),r(ldr,*),u(*),v(*),wk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c                                                                       
c     Inout  Q       Real(ldq,n)      orthogonal matrix from QR         
c     In     ldq     Integer          leading dimension of Q            
c     In     n       Integer          order of Q and R                  
c     Inout  R       Real(ldr,n)      upper triangular matrix R from QR 
c     In     ldr     Integer          leading dimension of R            
c     In     u       Real(c)          vector u of size n                
c     In     v       Real(c)          vector v of size n                
c     Out    wk      Real(c)          workspace of size n               
c                                                                      
c     on return                                                         
c                                                                       
c        Q       Q is the matrix with orthonormal columns in a QR       
c                decomposition of the matrix B = A + u*v'               
c                                                                       
c        R       R is the upper triangular matrix in a QR               
c                decomposition of the matrix B = A + u*v'               
c                                                                       
c     Description                                                       
c                                                                       
c     The matrices Q and R are a QR decomposition of a square matrix    
c     A = Q*R.                                                          
c     Given Q and R, qrupdt computes a QR decomposition of the rank one 
c     modification B = A + u*trans(v) of A. Here u and v are vectors.   
c                                                                       
c     Source : Algorithm 686                                            
c              Fortran subroutines for updating the QR decomposition    
c              ACM Transactions on Mathematical software, dec. 1990     
c                                                                       
c-----------------------------------------------------------------------------

c     Local variables and functions

      integer k,i
      double precision  c,s,tmp
      double precision  ddot

c     calculate wk = trans(Q)*u

      do i=1,n
         wk(i) = ddot(n,q(1,i),1,u,1)
      enddo

c     zero components wk(n),wk(n-1)...wk(2)
c     and apply rotators to R and Q.

      do k=n-1,1,-1
         call nuvgiv(wk(k),wk(k+1),c,s)
         call drot(n-k+1,r(k,k),ldr,r(k+1,k),ldr,c,s)
         call drot(n    ,q(1,k),1  ,q(1,k+1),1  ,c,s)
      enddo 

c     r(1,1:n) += wk(1)*v(1:n)      
      call daxpy(n,wk(1),v,1,r(1,1),ldr)
      
c     R is of upper hessenberg form. Triangularize R.

      do k=1,n-1
         call nuvgiv(r(k,k),r(k+1,k),c,s)
         call drot(n-k,r(k,k+1),ldr,r(k+1,k+1),ldr,c,s)
         call drot(n  ,q(1,k)  ,1  ,q(1,k+1)  ,1  ,c,s)
      enddo

      return
      end

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
