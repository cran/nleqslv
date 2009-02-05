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

      do 10 i=1,n
         wk(i) = ddot(n,q(1,i),1,u,1)
   10 continue

c     zero components wk(n),wk(n-1)...wk(2)
c     and apply rotators to R and Q.

      do 20 k=n-1,1,-1
         call nuvgiv(wk(k),wk(k+1),c,s)
         call drot(n-k+1,r(k,k),ldr,r(k+1,k),ldr,c,s)
         call drot(n    ,q(1,k),1  ,q(1,k+1),1  ,c,s)
   20 continue

c     r(1,1:n) += wk(1)*v(1:n)      
      call daxpy(n,wk(1),v,1,r(1,1),ldr)
      
c     R is of upper hessenberg form. Triangularize R.

      do 40 k=1,n-1
         call nuvgiv(r(k,k),r(k+1,k),c,s)
         call drot(n-k,r(k,k+1),ldr,r(k+1,k+1),ldr,c,s)
         call drot(n  ,q(1,k)  ,1  ,q(1,k+1)  ,1  ,c,s)
   40 continue

      return
      end
