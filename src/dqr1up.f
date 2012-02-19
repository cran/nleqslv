c The routines in this file come from the opensource library qrupdate
c (http://sourceforge.net/projects/qrupdate/).
c From that library the files dqr1up.f, dqrqh.f, dqhqr.f, dqrot.f, dqrtv1.f and dch1up.f
c have been combined in this file. No further changes have been made.
c 2012 Berend H. Hasselman

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dqr1up(m,n,k,Q,ldq,R,ldr,u,v,w)
c purpose:      updates a QR factorization after rank-1 modification
c               i.e., given a m-by-k orthogonal Q and m-by-n upper
c               trapezoidal R, an m-vector u and n-vector v,
c               this subroutine updates Q -> Q1 and R -> R1 so that
c               Q1*R1 = Q*R + u*v', and Q1 is again orthonormal
c               and R1 upper trapezoidal.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix Q.
c n (in)        number of columns of the matrix R.
c k (in)        number of columns of Q, and rows of R. Must be
c               either k = m (full Q) or k = n < m (economical form).
c Q (io)        on entry, the orthogonal m-by-k matrix Q.
c               on exit, the updated matrix Q1.
c ldq (in)      the leading dimension of Q. ldq >= m.
c R (io)        on entry, the upper trapezoidal m-by-n matrix R..
c               on exit, the updated matrix R1.
c ldr (in)      the leading dimension of R. ldr >= k.
c u (io)        the left m-vector. On exit, if k < m, u is destroyed.
c v (io)        the right n-vector. On exit, v is destroyed.
c w (out)       a workspace vector of size 2*k
c
      integer m,n,k,ldq,ldr
      double precision Q(ldq,*),R(ldr,*),u(*),v(*),w(*)
      external dqrqh,dqhqr,dqrot,dqrtv1
      external daxpy,ddot,dnrm2,dlamch,dscal,drot
      double precision ddot,dnrm2,dlamch,ru,ruu
      integer info,i
      logical full
c quick return if possible.
      if (k == 0 .or. n == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (k /= m .and. (k /= n .or. n > m)) then
        info = 3
      else if (ldq < m) then
        info = 5
      else if (ldr < k) then
        info = 7
      endif
      if (info /= 0) then
        call xerbla('DQR1UP',info)
        return
      end if

      full = k == m
c in the non-full case, we shall need the norm of u.
      if (.not.full) ru = dnrm2(m,u,1)
c form Q'*u. In the non-full case, form also u - Q*Q'u.
      do i = 1,k
        w(i) = ddot(m,Q(1,i),1,u,1)
        if (.not.full) call daxpy(m,-w(i),Q(1,i),1,u,1)
      end do
c generate rotations to eliminate Q'*u.
      call dqrtv1(k,w,w(k+1))
c apply rotations to R.
      call dqrqh(k,n,R,ldr,w(k+1),w(2))
c apply rotations to Q.
      call dqrot('B',m,k,Q,ldq,w(k+1),w(2))
c update the first row of R.
      call daxpy(n,w(1),v,1,R(1,1),ldr)
c retriangularize R.
      call dqhqr(k,n,R,ldr,w(k+1),w)
c apply rotations to Q.
      call dqrot('F',m,min(k,n+1),Q,ldq,w(k+1),w)
c in the full case, we're finished
      if (full) return
c compute relative residual norm
      ruu = dnrm2(m,u,1)
      ru = ru * dlamch('e')
      if (ruu <= ru) return
c update the orthogonal basis.
      call dscal(n,ruu,v,1)
      call dscal(m,1d0/ruu,u,1)
      call dch1up(n,R,ldr,v,w(k+1))
      do i = 1,n
        call drot(m,Q(1,i),1,u,1,w(k+i),v(i))
      end do
      end subroutine

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dqrqh(m,n,R,ldr,c,s)
c purpose:      brings an upper trapezoidal matrix R into upper
c               Hessenberg form using min(m-1,n) Givens rotations.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix R
c n (in)        number of columns of the matrix R
c R (io)        on entry, the upper Hessenberg matrix R
c               on exit, the updated upper trapezoidal matrix
c ldr (in)      leading dimension of R, >= m
c c(in)         rotation cosines, size at least min(m-1,n)
c s(in)         rotation sines, size at least min(m-1,n)
c
      integer m,n,ldr
      double precision R(ldr,*),c(*),s(*)
      external xerbla
      double precision t
      integer info,i,ii,j
c quick return if possible.
      if (m == 0 .or. m == 1 .or. n == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (ldr < m) then
        info = 4
      end if
      if (info /= 0) then
        call xerbla('DQRQH',info)
        return
      end if
      do i = 1,n
        ii = min(m-1,i)
c apply stored rotations, column-wise
        t = R(ii+1,i)
        do j = ii,1,-1
          R(j+1,i) = c(j)*t - s(j)*R(j,i)
          t = c(j)*R(j,i) + s(j)*t
        end do
        R(1,i) = t
      end do
      end subroutine

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dqhqr(m,n,R,ldr,c,s)
c purpose:      given an m-by-n upper Hessenberg matrix R, this
c               subroutine updates R to upper trapezoidal form
c               using min(m-1,n) Givens rotations.
c               (real version)
c arguments:
c m (in)        number of rows of the matrix R
c n (in)        number of columns of the matrix R
c R (io)        on entry, the upper Hessenberg matrix R
c               on exit, the updated upper trapezoidal matrix
c ldr (in)      leading dimension of R, >= m
c c(out)        rotation cosines, size at least min(m-1,n)
c s(out)        rotation sines, size at least min(m-1,n)
c
      integer m,n,ldr
      double precision R(ldr,*),c(*),s(*)
      external xerbla,dlartg
      double precision t
      integer info,i,ii,j
c quick return if possible.
      if (m == 0 .or. m == 1 .or. n == 0) return
c check arguments.
      info = 0
      if (m < 0) then
        info = 1
      else if (n < 0) then
        info = 2
      else if (ldr < m) then
        info = 4
      end if
      if (info /= 0) then
        call xerbla('DQHQR',info)
        return
      end if
      do i = 1,n
c apply stored rotations, column-wise
        t = R(1,i)
        ii = min(m,i)
        do j = 1,ii-1
          R(j,i) = c(j)*t + s(j)*R(j+1,i)
          t = c(j)*R(j+1,i) - s(j)*t
        end do
        if (ii < m) then
c generate next rotation
          call dlartg(t,R(ii+1,i),c(i),s(i),R(ii,i))
          R(ii+1,i) = 0d0
        else
          R(ii,i) = t
        end if
      end do
      end subroutine

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dqrot(dir,m,n,Q,ldq,c,s)
c purpose:      Apply a sequence of inv. rotations from right
c
c arguments:
c dir (in)      if 'B' or 'b', rotations are applied from backwards
c               if 'F' or 'f', from forwards.
c m (in)        number of rows of matrix Q
c n (in)        number of columns of the matrix Q
c Q (io)        on entry, the matrix Q
c               on exit, the updated matrix Q1
c ldq (in)      the leading dimension of Q
c c (in)        n-1 rotation cosines
c s (in)        n-1 rotation sines
c
      character dir
      integer m,n,ldq
      double precision Q(ldq,*),c(*),s(*)
      external drot,lsame
      logical lsame,fwd
      integer info,i
c quick return if possible
      if (m == 0 .or. n == 0 .or. n == 1) return
c check arguments.
      info = 0
      fwd = lsame(dir,'F')
      if (.not.(fwd .or. lsame(dir,'B'))) then
        info = 1
      else if (m < 0) then
        info = 2
      else if (n < 0) then
        info = 3
      else if (ldq < m) then
        info = 5
      end if
      if (info /= 0) then
        call xerbla('DQROT',info)
        return
      end if

      if (fwd) then
        do i = 1,n-1
          call drot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
      else
        do i = n-1,1,-1
          call drot(m,Q(1,i),1,Q(1,i+1),1,c(i),s(i))
        end do
      end if
      end subroutine

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dqrtv1(n,u,w)
c purpose:      generates a sequence of n-1 Givens rotations that
c               eliminate all but the first element of a vector u.
c arguments:
c n (in)        the length of the vector u
c u (io)        on entry, the vector u.
c               on exit, u(2:n) contains the rotation sines, u(1)
c               contains the remaining element.
c w (o)         on exit, w contains the rotation cosines.
c
      integer n
      double precision u(*),w(*)
      external dlartg
      double precision rr,t
      integer i
c quick return if possible.
      if (n <= 0) return
      rr = u(n)
      do i = n-1,1,-1
        call dlartg(u(i),rr,w(i),u(i+1),t)
        rr = t
      end do
      u(1) = rr
      end subroutine

c Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic
c
c Author: Jaroslav Hajek <highegg@gmail.com>
c
c This file is part of qrupdate.
c
c qrupdate is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 3 of the License, or
c (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this software; see the file COPYING.  If not, see
c <http://www.gnu.org/licenses/>.
c
      subroutine dch1up(n,R,ldr,u,w)
c purpose:      given an upper triangular matrix R that is a Cholesky
c               factor of a symmetric positive definite matrix A, i.e.
c               A = R'*R, this subroutine updates R -> R1 so that
c               R1'*R1 = A + u*u'
c               (real version)
c arguments:
c n (in)        the order of matrix R
c R (io)        on entry, the upper triangular matrix R
c               on exit, the updated matrix R1
c ldr (in)      leading dimension of R. ldr >= n.
c u (io)        the vector determining the rank-1 update
c               on exit, u contains the rotation sines
c               used to transform R to R1.
c w (out)       cosine parts of rotations.
c
      integer n,ldr
      double precision R(ldr,*),u(*)
      double precision w(*)
      external dlartg
      double precision rr,ui,t
      integer i,j

      do i = 1,n
c apply stored rotations, column-wise
        ui = u(i)
        do j = 1,i-1
          t = w(j)*R(j,i) + u(j)*ui
          ui = w(j)*ui - u(j)*R(j,i)
          R(j,i) = t
        end do
c generate next rotation
        call dlartg(R(i,i),ui,w(i),u(i),rr)
        R(i,i) = rr
      end do
      end subroutine
