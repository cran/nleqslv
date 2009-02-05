
      subroutine liqrfa(a, lda, n, tau, work, info)
      integer  lda, n, info
      double precision  a(lda,*), tau(*), work(*)

c-------------------------------------------------------------
c
c     QR decomposition of A(n,n)
c
c     Arguments
c
c      Inout A        Real(Lda,n)    Matrix to transform.
c
c      In    lda      Integer        Leading dimension of A
c      In    n        Integer        number of rows/cols A
c
c      Out   tau      Real(n)        Information for recovering
c                                    orthogonal part of decomposition
c      Out   work     Real(*)        workspace     
c
c     Lapack unblocked QR
c
c-------------------------------------------------------------

      call dgeqr2(n,n,A,lda,tau,work,info)

      return
      end

c=============================================================

      subroutine liqrqt(a, lda, n, tau, qty, work, info)
      integer lda, n, info
      double precision a(lda,*), tau(*), qty(*), work(*)

c-------------------------------------------------------------
c      Arguments
c
c      In    A     Real(Lda, n)    QR decomposition
c
c      In    Lda   Integer         Leading dimension A
c      In    n     Integer         Number of rows/columns in A
c
c      In    tau   Integer         Householder constants from QR
c
c      Inout qty   Real(n)         On input, vector y
c                                  On output, trans(Q)*y
c
c      Out   work  Real(n)         workspace
c                                   dimension as called here with side='L'
c
c     Liqrqt calculates trans(Q)*y from the QR decomposition
c
c     Lapack unblocked
c-------------------------------------------------------------

      call dorm2r('L','T',n,1,n,A,lda,tau,qty,n,work,info)

      return
      end

c=============================================================

      subroutine liqrqq(q,ldq,tau,n,work,info)
      integer n,ldq, info
      double precision  q(ldq,*),tau(*),work(*)

c     Arguments
c
c     Inout    Q     Real(ldq,*)     On input, QR decomposition
c                                    On output, the full Q
c     In       ldq   Integer         leading dimension of Q
c     In       tau   Real(n)         Householder constants of
c                                     the QR decomposition
c     In       n     Integer         number of rows/columns in Q
c     Out      work  Real(n)         workspace of length n
c     
c     Generate Q from QR decomposition Liqrfa (dgeqr2)
c
c     Lapack unblocked
c-------------------------------------------------------------
 
      call dorg2r(n,n,n,q,ldq,tau,work,info)    
      
      return
      end
      
