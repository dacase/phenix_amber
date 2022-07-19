#include "../include/dprec.fh"

      subroutine wallclock( wallc )
      implicit none
      _REAL_, intent(out) ::  wallc
      integer, save :: ncalls = 0
      integer :: count, rate, n

      call system_clock( COUNT=count, COUNT_RATE=rate)
      wallc = dble(count)/dble(rate)
      ncalls = ncalls + 1
      return
 
      entry nwallclock ( n )
      n = ncalls
      end
