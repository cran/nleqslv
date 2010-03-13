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

      subroutine liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
      integer n,ldr
      double precision  r(ldr,*),diag(*),b(*),x(*),sdiag(*),wrk(*)

c-----------------------------------------------------------------------------
c
c     Arguments
c
c       In     n      Integer         order of R.
c       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
c                                     unaltered
c                                     strict lower triangle contains 
c                                        transposed strict upper triangle of the upper
c                                        triangular matrix S.
c
c       In     diag   Real(*)         vector with diagonal elements of matrix D
c                                     all elements > 0
c
c       In     ldr    Integer         leading dimension of the array R.
c       In     b      Real(*)         vector of size n
c
c       Out    x      Real(*)         vector of size n
c                                     on output contains least squares solution
c                                     of the system R*x = b, D*x = 0.
c
c       Out    sdiag  Real(*)         vector of size n, containing the
c                                     diagonal elements of the upper
c                                     triangular matrix S.
c
c       Out    wrk    Real(*)         workspace of size n.
c
c     Description
c
c     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
c     and an n-vector b, determine an x which solves the system
c
c         |R*x| = |b|
c         |D*x| = |0|
c
c     in the least squares sense.
c     The routine is used when the matrix R from the QR decomposition of a Jacobian
c     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
c     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
c     diagonal matrix. The idea is then to solve (J^T * J + mu*I)x=b where mu
c     is a small positive number. Using a QR decomposition of J solving this system
c     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
c     Solving this system is equivalent to solving the above least squares problem.
c     On output the routine also provides an upper triangular matrix S such that
c     (see description of arguments above for the details)
c
c         (trans(R)*R + D*D) = trans(S)*S .
c
c     Method used here is described in
c     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
c     page 259--261
c-----------------------------------------------------------------------------

      integer j,k
      double precision  bj,c,s,sum,temp
      double precision  ddot
      double precision Rzero, Rone
      parameter(Rzero=0.0d0, Rone=1.0d0)

c     copy R and b to preserve input and initialise S.
c     Save the diagonal elements of R in wrk.
c     Beware: the algorithm operates on an upper triangular matrix,
c     which is stored in lower triangle of R.
c     
      do j=1,n
         call dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
         wrk(j) = r(j,j)
      enddo
      call dcopy(n,b,1,x,1)

c     eliminate the diagonal matrix D using givens rotations.
c     Nocedal method: start at the bottom right
c     after 100 loop had finished R contains diagonal of S
c     save in sdiag and restore original diagonal of R

      do j=n,1,-1

c        initialise the row of D to be eliminated

         call nuzero(n-j+1,sdiag(j))
         sdiag(j) = diag(j)

c        the transformations to eliminate the row of D

         bj = Rzero
         do k=j,n

c           determine a givens rotation which eliminates the
c           appropriate element in the current row of D.
c           accumulate the transformation in the row of S.

c           eliminate the diagonal element in row j of D
c           this generates fill-in in columns [j+1 .. n] of row j of D
c           successively eliminate the fill-in with givens rotations
c           for R[j+1,j+1] and D[j,j+1].
c           rows of R have been copied into the columns of R in 10 loop
c           perform all operations on those columns to preserve the original R

            if (sdiag(k) .ne. Rzero) then

               call nuvgiv(r(k,k),sdiag(k),c,s)
               if( k .lt. n ) then
                   call drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
               endif

c              compute the modified element of (b,0).

               temp =  c*x(k) + s*bj
               bj   = -s*x(k) + c*bj
               x(k) = temp

            endif

         enddo

      enddo

c     retrieve diagonal of S from diagonal of R
c     restore original diagonal of R

      do k=1,n
         sdiag(k) = r(k,k)
         r(k,k) = wrk(k)
      enddo

c     x now contains modified b
c     solve trans(S)*x = x

      x(n) = x(n) / sdiag(n)
      do j=n-1,1,-1
         sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
         x(j) = (x(j) - sum)/sdiag(j)
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
