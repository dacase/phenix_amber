
  else if ( eedmeth == 2 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 2
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do

    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1,icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Linear lookup on switch:
      delrinv = one / sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      xx = eedtbdns*x + 1
      ind = int(xx)
      dx = xx - ind
      switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
      d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
        if (ifcr .ne. 0) then
          b2 = b0 * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      comm1 = cgi*cgj
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End prologue loop

    if (tvips) then

      ! Use IPS for L-J energy:
#     include "ips_lj.h"
    else

      ! regular epilogue:
#     include "ew_directe.h"
    end if

    ! Now loop over the 12-10 LJ terms for eedmeth = 2
    icount = 0
    do m = nvdw + 1, ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Linear lookup on switch:
      delrinv = one/sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      xx = eedtbdns*x + 1
      ind = int(xx)
      dx = xx - ind
      switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
      d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi*cgj
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
        if (ifcr .ne. 0) then
          b2 = lfac * b0
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      comm1 = cgi*cgj
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor( i, b0*cgj )
        call cr_add_dcdr_factor( j, b0*cgi )
      end if
#endif
      ecur = comm1*b0        
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End epilogue loop

#   include "ew_directe2.h"

  else if ( eedmeth == 3 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 3
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do
      
    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1,icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Explicit function call:
      delrinv = one/sqrt(delr2)
      delr = delr2*delrinv
      delr2inv = delrinv*delrinv
      x = dxdr * delr
      call get_ee_func(x, switch, d_switch_dx, ee_type)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
        if (ifcr .ne. 0) then
          b2 = b0 * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      comm1 = cgi * cgj
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new)=delr2inv
      cache_df(im_new)=dfee
    end do
    ! End electrostatic loop

    if (tvips) then

      ! Use IPS for L-J energy:
#     include "ips_lj.h"
    else

      ! regular epilogue:
#     include "ew_directe.h"
    end if

    ! Now loop over the 12-10 LJ terms for eedmeth = 3
    icount = 0
    do m = nvdw + 1, ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Explicit function call:
      delrinv = one / sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      call get_ee_func(x, switch, d_switch_dx, ee_type)
      b0 = switch*delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi*cgj*lfac
        if (ifcr .ne. 0) then
          b2 = b0 * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      comm1 = cgi * cgj
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End epilogue loop

#   include "ew_directe2.h"
  else if (eedmeth == 4) then

    ! Loop over the 12-6 LJ terms for eedmeth = 4
    icount = 0
    do m = 1, nvdw
#     include "ew_directp.h"
    end do
      
    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, delrinv*cgj)
          call cr_add_dcdr_factor(j, delrinv*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi*cgj*lfac*delrinv
        if (ifcr .ne. 0) then
          b2 = lfac * delrinv
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      b0 = cgi * cgj * delrinv
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, delrinv*cgj)
        call cr_add_dcdr_factor(j, delrinv*cgi)
      end if
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = b0 * delr2inv
         
      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End prologue loop

    if (tvips) then

      ! Use IPS for L-J energy:
#     include "ips_lj.h"
    else

      ! regular epilogue:
#     include "ew_directe.h"
    end if

    ! Now loop over the 12-10 LJ terms for eedmeth = 4
    icount = 0
    do m = nvdw+1, nvdw+nhbnd
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount

      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, delrinv*cgj)
          call cr_add_dcdr_factor(j, delrinv*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delrinv
        if (ifcr .ne. 0) then
          b2 = lfac * delrinv
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      b0 = cgi * cgj * delrinv
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, delrinv*cgj)
        call cr_add_dcdr_factor(j, delrinv*cgi)
      end if
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = b0 * delr2inv
      cache_r2(im_new)=delr2inv
      cache_df(im_new)=dfee
    end do

#   include "ew_directe2.h"

#ifdef MPI /* SOFT CORE */
    ! Now loop over the softcore (modified 12-6) LJ terms for eedmeth = 4
    icount = 0
    do m = nvdw+nhbnd+1,ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#  if defined(LES) 
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, delrinv*cgj)
          call cr_add_dcdr_factor(j, delrinv*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delrinv
        if (ifcr .ne. 0) then
          b2 = lfac * delrinv
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#  else
      b0 = cgi*cgj*delrinv
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, delrinv*cgj)
        call cr_add_dcdr_factor(j, delrinv*cgi)
      end if
#  endif /* LES */
      eelt = eelt + b0
      dfee = b0 * delr2inv
      cache_r2(im_new) = delr2

      ! Contrary to the 12-6 and 12-10 cases above, cache_r2 contains r^2 here
      cache_df(im_new) = dfee
    end do

    ! V1 uses ew_directe3.h, in which softcore atoms are treated as appearing,
    ! i.e. fully interacting at lambda=1 and 'soft' at small lambda
    ! V0 uses ew_directe4.h, in which softcore atoms are treated as vanishing, 
    ! i.e. fully interacting at lambda=0 and 'soft' at large lambda
    if (isProcessV1) then
#     include "ew_directe3.h"
    else
#     include "ew_directe4.h"
    end if
#endif /* MPI for SOFT CORE */

  else if ( eedmeth == 5 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 5
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do
      
    !  Calculation starts: loop over the data gathered in the temporary
    !  array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Use a distance-dependent dielectric of 1/r:
      delr2inv = one / delr2
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delr2inv
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, delr2inv*cgj)
          call cr_add_dcdr_factor(j, delr2inv*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delr2inv
        if (ifcr .ne. 0) then
          b2 = delr2inv * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      b0 = cgi*cgj*delr2inv
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, delr2inv*cgj)
        call cr_add_dcdr_factor(j, delr2inv*cgi)
      end if
#endif 
      ecur = b0
      eelt = eelt + ecur
      dfee = two * b0 * delr2inv

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    if (tvips) then

      ! Use IPS for L-J energy:
#     include "ips_lj.h"
    else

      ! regular epilogue:
#     include "ew_directe.h"
    end if

    !     Now loop over the 12-10 LJ terms for eedmeth = 5
    icount = 0
    do m = nvdw+1, ntot
#     include "ew_directp.h"
    end do

    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Use dielectric of 1/r:
      delr2inv = one / delr2
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delr2inv
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, delr2inv*cgj)
          call cr_add_dcdr_factor(j, delr2inv*cgi)
        end if
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delr2inv
        if (ifcr .ne. 0) then
          b2 = delr2inv * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
      end if
#else
      b0 = cgi * cgj * delr2inv
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, delr2inv*cgj)
        call cr_add_dcdr_factor(j, delr2inv*cgi)
      end if
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = two * b0 * delr2inv
      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee

    end do
#   include "ew_directe2.h"
  else if (eedmeth == 6) then

    ! Loop over the 12-6 LJ terms for eedmeth = 6
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do
    call vdinvsqrt( icount, cache_r2, cache_df )

    ! the df cache now stores the reciprocal of delta r, variable delrinv.
    if (teips .or. tvips) then
      nnbips = nnbips + icount*2
    end if
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
      delrinv = cache_df(im_new)
      delr2inv = delrinv * delrinv

      ! Use the ips potential:
      uips = ripsr * delr2 * delrinv
      uips2 = rips2r * delr2
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
        b1 = delrinv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b1 = lfac * delrinv
        b0 = cgi * cgj * b1
      end if
#else
      b0 = cgi * cgj * ripsr
      b1 = ripsr
#endif
      pipse = one/uips + AIPSE(0) + &
              uips2*(AIPSE(1) + uips2*(AIPSE(2) + uips2*AIPSE(3)))
      dpipse = -one/uips + &
               uips2*(BIPSE(1) + uips2*(BIPSE(2) + uips2*BIPSE(3)))
      ecur = b0*(pipse - pipsec)
      eelt = eelt + ecur
      dfee = -b0 * dpipse * delr2inv
      if (ifcr .ne. 0) then
        b2 = b1*(pipse - pipsec)
        call cr_add_dcdr_factor(i, b2*cgj)
        call cr_add_dcdr_factor(j, b2*cgi)
      end if
      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do

    if (TVIPS) then

      ! Use IPS for L-J energy
#     include "ips_lj.h"
    else

      ! epilogue
#     include "ew_directe.h"
    endif

    ! Now loop over the 12-10 LJ terms for eedmeth = 6
    icount = 0
    do m = nvdw+1, ntot
#     include "ew_directp.h"
    end do
    call vdinvsqrt(icount, cache_r2, cache_df)
    if (teips .or. tvips) then
      nnbips = nnbips + icount*2
    end if
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
      delrinv = cache_df(im_new)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi*cgj*delrinv
        b1 = delrinv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b1 = lfac * delrinv
        b0 = cgi * cgj * b1
      end if
#else
      b0 = cgi * cgj * ripsr
      b1 = ripsr
#endif

      ! Use a ips potential:
      uips = ripsr * delr2 * delrinv
      uips2 = rips2r * delr2
      pipse = one/uips + AIPSE(0) + &
              uips2*(AIPSE(1) + uips2*(AIPSE(2) + uips2*AIPSE(3)))
      dpipse = -one/uips + &
               uips2*(BIPSE(1) + uips2*(BIPSE(2) + uips2*BIPSE(3)))
      ecur = b0*(pipse - pipsec)
      eelt = eelt + ecur
      if (ifcr /= 0) then
        b2 = b1*(pipse - pipsec)
      call cr_add_dcdr_factor(i, b2*cgj)
      call cr_add_dcdr_factor(j, b2*cgi)
    end if
    dfee = -b0 * dpipse * delr2inv

    cache_df(im_new)=dfee
  end do

#   include "ew_directe2.h"
