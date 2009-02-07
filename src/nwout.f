
c-----------------------------------------------------------------------

      subroutine nwstrot(s)
      character*(*) s

      integer slen
      slen = len_trim(s)
      call nwstrot0(s,slen)

      return
      end

      subroutine nwckot(i,j,aij,wi)
      integer i,j
      double precision aij,wi

      character*80 s

c     error message for check analytic jacobian

      write(s,900) i,j
      call nwstrot(s)
      write(s,901) aij,wi
      call nwstrot(s)

 900  format('Chkjac  possible error in jacobian(',i4,',',i4,')')
 901  format('       ', d20.13, ' Estimated = ', d20.13)

      return
      end
