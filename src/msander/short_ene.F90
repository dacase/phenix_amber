!<compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

#if 1  /* always use the non-OpenMP code for now */
!  #if defined(LES) || defined(MPI) /* use the older, non-openMP version */

!------------------------------------------------------------------------------
! get_nb_energy: the main routine for vdw, hbond, and direct space Ewald sum
!                computations.  It is structured for parallelism.
!
! Arguments:
!
!------------------------------------------------------------------------------
subroutine get_nb_energy(iac, ico, ntypes, charge, cn1, cn2, cn6, force, &
                         numatoms, ipairs, ewaldcof, eedtbdns, eed_cub, &
                         eed_lin, maxnblst, eelt, evdw, ehb, dir_vir, eedvir, &
                         filter_cut, ee_type, eedmeth, dxdr, cn3, &
                         cn4, cn5, epol, dipole, field)

  use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
                     myindexlo, myindexhi, numimg
  use constants, only : zero
  use stack
#ifdef MPI /* SOFT CORE */
  use softcore, only : sc_ener
  use nblist, only : numsc
#endif

  implicit none
  character(kind=1, len=13) :: routine="get_nb_energy"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
#include "flocntrl.h"
#include "def_time.h"
   
  integer l_real_df, l_real_x, l_real_y, l_real_z, l_real_r2, l_int
  integer numatoms, maxnblst
  integer iac(*), ico(*), ntypes, ee_type, eedmeth
  _REAL_ charge(*), cn1(*), cn2(*), cn6(*)
  _REAL_ ewaldcof, eedtbdns, dxdr, eed_cub(4,*), eed_lin(2,*), dir_vir(3,3)
  integer ipairs(maxnblst)
  _REAL_ force(3,numatoms), eelt, epol, evdw, ehb
  _REAL_ eedvir, filter_cut, dipole(3,*), field(3,*)
  _REAL_ cn3(*), cn4(*), cn5(*)
   
  integer index, numpack, i, k, ncell_lo, ncell_hi, ntot, nvdw, nhbnd
  _REAL_ xk(3)
  integer ncache
   
  if (do_dir == 0) then
    return
  end if
  eelt = zero
  epol = zero
  evdw = zero
  ehb = zero
  eedvir = zero
#ifdef MPI /* SOFT CORE */
  sc_ener(7) = 0.0d0
  sc_ener(8) = 0.0d0  
  sc_ener(9) = 0.0d0
  sc_ener(10) = 0.d0
  sc_ener(11) = 0.d0
#endif
  dir_vir(1:3,1:3) = zero
  numpack = 1
  call timer_start(TIME_SHORT_ENE)
  do index = myindexlo, myindexhi
    if (numimg(index) > 0) then
      ncell_lo = nlogrid(index)
      ncell_hi = nhigrid(index)
      do k = ncell_lo, ncell_hi
        i = bckptr(k)
#ifdef MPI /* SOFT CORE */
        ! SOFT CORE contribution in numsc
        ntot = numvdw(i) + numhbnd(i) + numsc(i)
#else
        ntot = numvdw(i) + numhbnd(i)
#endif
        nvdw = numvdw(i)
        nhbnd = numhbnd(i)
        if (ntot > 0) then
          xk(:) = imagcrds(:,k)

          ! Allocate 6 temporary caches for performance optimizations
#ifdef MPI /* SOFT CORE */
          if (nhbnd > numsc(i)) then
#endif
            ncache = max(nvdw, nhbnd)
#ifdef MPI /* SOFT CORE */
          else
            ncache = max( nvdw, numsc(i) )
          end if
#endif
          call get_stack(l_real_df, ncache, routine)
          call get_stack(l_real_x, ncache, routine)
          call get_stack(l_real_y, ncache, routine)
          call get_stack(l_real_z, ncache, routine)
          call get_stack(l_real_r2, ncache, routine)
          call get_istack(l_int, ncache, routine)
          if (.not. rstack_ok) then
            deallocate(r_stack)
            allocate(r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
          endif
          if (.not. istack_ok) then
            deallocate(i_stack)
            allocate(i_stack(1:lastist), stat=alloc_ier)
            call reassign_istack(routine)
          endif
          REQUIRE(rstack_ok)
          REQUIRE(istack_ok)
          call short_ene(i, xk, ipairs(numpack), ntot, nvdw, nhbnd, &
                         eedtbdns, eed_cub, eed_lin, charge, ntypes, iac, &
                         ico, cn1, cn2, cn6, filter_cut, eelt, evdw, force, &
                         dir_vir, ee_type, eedmeth, dxdr, eedvir, &
                         r_stack(l_real_df), r_stack(l_real_x), &
                         r_stack(l_real_y), r_stack(l_real_z), &
                         r_stack(l_real_r2), i_stack(l_int), cn3, cn4, cn5)
                      
          call free_stack(l_real_r2,routine)
          call free_stack(l_real_z,routine)
          call free_stack(l_real_y,routine)
          call free_stack(l_real_x,routine)
          call free_stack(l_real_df,routine)
          call free_istack(l_int,routine)
        end if  ! ( ntot > 0 )
        numpack = numpack + ntot
      end do  !  k = ncell_lo,ncell_hi
    end if  ! ( numimg(k) > 0 )
  end do
  ! nd loop over this process's assigned atoms

  call timer_stop(TIME_SHORT_ENE)

  return

end subroutine get_nb_energy 

!------------------------------------------------------------------------------
! short_ene: calculate the direct space component of the Ewald potentials.
!            This subroutine was modified by Nathalie Godbout, sgi 04/26/00.
!            Additional optimizations aimed at IA32 SSE2 were provided by Scott
!            Brozell, TSRI Oct 2002.  More tweaks were made by Ross Walker,
!            TSRI 2005.
!------------------------------------------------------------------------------
subroutine short_ene(i, xk, ipairs, ntot, nvdw, nhbnd, eedtbdns, &
                     eed_cub, eed_lin, charge, ntypes, iac, ico, cn1, cn2, &
                     cn6, filter_cut, eelt, evdw, force, dir_vir, ee_type, &
                     eedmeth, dxdr, eedvir, cache_df, cache_x, cache_y, &
                     cache_z, cache_r2, cache_bckptr, cn3, cn4, cn5)
  use nblist, only: imagcrds, bckptr, tranvec, cutoffnb, volume
  use constants, only: zero, one, two, half, third, TWOPI, six, twelve
  use file_io_dat
#ifdef LES
  use les_data, only: cnum, lestyp, lestmp, lesfac, lfac, nlesty
#endif
  use nbips, only: teips, tvips, nnbips, rips2, ripsr, rips2r, rips6r, &
                   rips12r, aipse, aipsvc, aipsva, bipse, bipsvc, bipsva, &
                   pipsec, pipsvcc, pipsvac
#ifdef MPI /* SOFT CORE */
  use softcore, only: scalpha, scbeta, sigma6, foureps, sc_dvdl, &
                      isProcessV1, sc_ener, nsc, oneweight, weight0, &
                      weight1, sceeorder, sc_dvdl_ee
  use mbar, only: ifmbar, bar_i, bar_states, bar_lambda, bar_cont, do_mbar
#endif
  use crg_reloc, only: ifcr, cr_add_dcdr_factor

  implicit none
  _REAL_ xk(3)
  integer i, nvdw, nhbnd, ntot
  integer ipairs(*), ee_type, eedmeth
  _REAL_ eed_cub(*), eed_lin(2,*), charge(*), dir_vir(3,3), eedvir
  _REAL_ eedtbdns, filter_cut, dxdr
  integer ntypes, iac(*), ico(*)
  _REAL_ cn1(*), cn2(*), cn6(*), eelt, evdw, force(3,*)
  integer ic,j,m,n,ind,iaci
  _REAL_ del, delrinv, delr12inv
  _REAL_ switch, d_switch_dx
  _REAL_ b0, b1, b2
  _REAL_ ee_vir_iso
  _REAL_ filter_cut2, xx
  _REAL_ comm1
  _REAL_ xktran(3,18)
  _REAL_ e3dx, e4dx
  _REAL_ cn3(*), cn4(*), cn5(*)
  _REAL_ ecur
  integer, parameter :: mask27 = 2**27 - 1
#ifdef LES
#  include "ew_cntrl.h"
#endif
  _REAL_ delx,dely,delz,delr, delr2, cgi, cgj, delr2inv, r6, f6, f12, df, &
         dfee, dx, x, dfx, dfy, dumy, dfz, dumz, dumx
#ifdef MPI
  _REAL_ denom, denom_n, delr_n, switch_c, denom2, denom3, rfour
#endif

#include "../include/md.h"
  _REAL_ uips, uips2, uips4, uips2r, uips6r, uips12r
  _REAL_ pipse, dpipse, pvc, dvcu, pva, dvau
  integer itran

  ! variables for conditionally cached data.
  _REAL_ cache_df(*), cache_x(*), cache_y(*)
  _REAL_ cache_z(*), cache_r2(*)
  integer im_new, icount
  integer cache_bckptr(*)
  _REAL_ time0, time1

  ee_vir_iso = zero
  del = one / eedtbdns
  dumx = zero
  dumy = zero
  dumz = zero
  filter_cut2 = filter_cut * filter_cut
  cgi = charge(i)
  iaci = ntypes * (iac(i) - 1)
   
  xktran(1,:) = tranvec(1,:) - xk(1)
  xktran(2,:) = tranvec(2,:) - xk(2)
  xktran(3,:) = tranvec(3,:) - xk(3)

#ifdef LES
  lestmp=nlesty*(lestyp(i)-1)
#endif
   
  ! The "eedmeth" decision is unrolled here, since this provides significant
  ! efficiency improvements, in spite of greatly increasing the code length.
  ! Each eedmeth case calculates the 12-6 LJ terms and 12-10 LJ terms in
  ! separate sets of fissioned loops.  Each loop set consists of three loop
  ! fission fragments: prologue, electrostatic, and epilogue.
  !
  ! The prologue loop computes delta r**2 and conditionally caches it:
  ! delta x, delta y, delta z, and the bckptr index.  The electrostatic loop
  ! calculates and caches the direct Ewald sum.  The epilogue loop computes
  ! the van der Waals interaction and updates the forces using the previously
  ! cached data.  This substantial code tuning is designed to enable software
  ! pipelining of loops for superscalar and similar architectures.
  !
  ! For eedmeth = 1 12-6 LJ, the most common case, the electrostatic loop is
  ! fissioned into two fragments. The first loop calculates and caches
  ! principally the reciprocal square root of delta r**2.  The second loop
  ! completes the electrostatic work.  This code tuning is designed to enable
  ! SIMD vectorization of the reciprocal square root on IA32 SSE2 compatible
  ! platforms.
  if (eedmeth == 1) then
     
    !-------------------------------------------------------------
    ! Loop over the 12-6 LJ terms for eedmeth = 1
    !-------------------------------------------------------------
    icount = 0
    do m = 1, nvdw
!   #     include "ew_directp.h"

    ! prologue: gather the data and put it in temporary arrays.

       n=ipairs(m)
       itran=ishft(n,-27)
       n = iand(n,mask27)
       j = bckptr(n)
       delx = imagcrds(1,n) + xktran(1,itran)
       dely = imagcrds(2,n) + xktran(2,itran)
       delz = imagcrds(3,n) + xktran(3,itran)
       delr2 = delx*delx + dely*dely+delz*delz
   
       if ( delr2 < filter_cut2 )then
          icount = icount + 1
          cache_x(icount) = delx
          cache_y(icount) = dely
          cache_z(icount) = delz
          cache_r2(icount) = delr2
          cache_bckptr(icount) = j
       end if
    end do

    !  calculation starts: loop over the data gathered in the temporary
    !  array caches.
    call vdinvsqrt(icount, cache_r2, cache_df)
    cache_r2(1:icount) = dxdr*cache_df(1:icount) * cache_r2(1:icount)

    ! The df cache now stores the reciprocal of delta r, variable delrinv.
    ! The r2 cache contains delta r times the derivative; this product,
    ! dxdr*delr, is denoted as variable x below.
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delrinv = cache_df(im_new)
      x = cache_r2(im_new)

      ! Cubic spline on switch:
      ind = int(eedtbdns*x)
      dx = x - ind*del
      ind = 4*ind
      e3dx = dx * eed_cub(3+ind)
      e4dx = dx * dx * eed_cub(4+ind)
      switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
               (e3dx + e4dx*third)*half)
      d_switch_dx = eed_cub(2+ind) + e3dx + e4dx*half
         
      ! Tom Darden got the idea for B_l from Walter Smith's CCP5 article 1982
      ! Ewald for point multipoles.
      b0 = switch*delrinv
      b1 = b0 - d_switch_dx*dxdr
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here
        comm1 = cgi * cgj

#if 0
        ! calculate contribution of dc/dr to force
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
#endif
      else
        lfac = lesfac(lestmp + lestyp(j))
        comm1 = cgi * cgj * lfac

#if 0
        ! Calculate contribution of dc/dr to force
        if (ifcr .ne. 0) then
          b2 = b0*lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
#endif
      end if
#else
      comm1 = cgi*cgj
 
#if 0
      ! Calculate contribution of dc/dr to force
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
#endif
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1*b1

      delr2inv = delrinv*delrinv
      dfee = dfee*delr2inv
      cache_r2(im_new)=delr2inv
      cache_df(im_new)=dfee
    end do
    ! End prologue loop

    if (tvips) then
      ! Use IPS for L-J energy:
#     include "ips_lj.h"
    else
        ! regular epilogue:
        !call wallclock(time0)
!   #     include "ew_directe.h"

! epilogue: 12-6 LF terms

        do im_new = 1,icount
           j = cache_bckptr(im_new)

           dfee = cache_df(im_new)
           delx = cache_x(im_new)
           dely = cache_y(im_new)
           delz = cache_z(im_new)
           delr2inv = cache_r2(im_new)

           ic = ico(iaci+iac(j))
           r6 = delr2inv*delr2inv*delr2inv
           delr12inv = r6 * r6
#ifdef LES 
           lfac=lesfac(lestmp+lestyp(j))
           f6 = cn2(ic)*r6*lfac
           f12 = cn1(ic)*delr12inv*lfac
#else
           f6 = cn2(ic)*r6
           f12 = cn1(ic)*delr12inv
#endif
           evdw = evdw + f12 - f6
           df = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv

           dfx = delx*df
           dfy = dely*df
           dfz = delz*df
           dumx = dumx + dfx
           dumy = dumy + dfy
           dumz = dumz + dfz
           force(1,j) = force(1,j) + dfx
           force(2,j) = force(2,j) + dfy
           force(3,j) = force(3,j) + dfz
        end do  !  im_new = 1,icount
        !call wallclock(time1)
        !time_ewd = time_ewd + time1-time0
    end if
      
    ! Now loop over the 12-10 LJ terms for eedmeth = 1
    icount = 0
    do m = nvdw+1,nvdw+nhbnd

! prologue: gather the data and put it in temporary arrays.

       n=ipairs(m)
       itran=ishft(n,-27)
       n = iand(n,mask27)
       j = bckptr(n)
       delx = imagcrds(1,n) + xktran(1,itran)
       dely = imagcrds(2,n) + xktran(2,itran)
       delz = imagcrds(3,n) + xktran(3,itran)
       delr2 = delx*delx + dely*dely+delz*delz
   
       if ( delr2 < filter_cut2 )then
          icount = icount + 1
          cache_x(icount) = delx
          cache_y(icount) = dely
          cache_z(icount) = delz
          cache_r2(icount) = delr2
          cache_bckptr(icount) = j
       end if
    end do

    call vdinvsqrt(icount, cache_r2, cache_df)
    cache_r2(1:icount) = dxdr * cache_df(1:icount) * cache_r2(1:icount)

    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delrinv = cache_df(im_new)
      x = cache_r2(im_new)
      
      ! Cubic spline on switch:
      ind = int(eedtbdns*x)
      dx = x - ind*del
      ind = 4*ind

      e3dx = dx*eed_cub(3+ind)
      e4dx = dx*dx*eed_cub(4+ind)
      switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + (e3dx + e4dx*third)*half)
      d_switch_dx = eed_cub(2+ind) + e3dx + e4dx*half
         
      ! Tom Darden Got the idea for B_l from Walter Smith's CCP5 article 1982
      ! Ewald for point multipoles
      b0 = switch * delrinv
      b1 = b0 - d_switch_dx*dxdr
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
            
        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here
        comm1 = cgi*cgj
#if 0
        if (ifcr .ne. 0) then
          call cr_add_dcdr_factor( i, b0*cgj )
          call cr_add_dcdr_factor( j, b0*cgi )
        end if
#endif
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
#if 0
        if (ifcr .ne. 0) then
          b2 = b0 * lfac
          call cr_add_dcdr_factor(i, b2*cgj)
          call cr_add_dcdr_factor(j, b2*cgi)
        end if
#endif
      end if
#else
      comm1 = cgi*cgj
#if 0
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
#endif
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      delr2inv = delrinv*delrinv
      dfee = dfee*delr2inv
      cache_df(im_new)=dfee
    end do
    ! End electrostatic loop

    ! epilogue: 12-10 LJ terms.

    !  --- this code allows 10-12 terms; in many (most?) (all?) cases, the
    !       only "nominal" 10-12 terms are on waters, where the asol and bsol
    !       parameters are always zero; hence we can skip this part.

    do im_new = 1,icount
       j = cache_bckptr(im_new)

       df =   cache_df(im_new)
       delx = cache_x(im_new)
       dely = cache_y(im_new)
       delz = cache_z(im_new)

       dfx = delx*df
       dfy = dely*df
       dfz = delz*df
       dumx = dumx + dfx
       dumy = dumy + dfy
       dumz = dumz + dfz
       force(1,j) = force(1,j) + dfx
       force(2,j) = force(2,j) + dfy
       force(3,j) = force(3,j) + dfz
    end do  !  im_new = 1,icount

!============================================================================
#ifdef MPI /* SOFT CORE */

    ! Now loop over the softcore (modified 12-6) LJ terms for eedmeth = 1
    icount = 0

    ! run over the 3rd subset of the pairlist
    do m = nvdw+nhbnd+1,ntot      
#     include "ew_directp.h"
    end do
      
    call vdinvsqrt(icount, cache_r2, cache_df)

    ! cache_r2 is saved for the vdW interactions below
    cache_r2(1:icount) = dxdr * cache_df(1:icount) * cache_r2(1:icount)

    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delrinv = cache_df(im_new)
      x = cache_r2(im_new)
         
      ! Cubic spline on switch:
      ind = int(eedtbdns*x)
      dx = x - ind*del
      ind = 4*ind
      e3dx = dx * eed_cub(3+ind)
      e4dx = dx * dx * eed_cub(4+ind)
      switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + (e3dx + e4dx*third)*half)
      d_switch_dx = eed_cub(2+ind) + e3dx + e4dx*half
         
      ! Tom Darden Got the idea for B_l from Walter Smith's CCP5 article 1982
      ! Ewald for point multipoles
      cgj = charge(j)
      comm1 = cgi * cgj
      delr2inv = delrinv * delrinv
      if (nsc(i) == nsc(j)) then          

        ! Both atoms are softcore and have normal Coulomb interactions.
        ! First remove the reciprocal contribution from the main energy array.
        switch_c = 1.0d0 - switch 
        eelt = eelt - comm1 * delrinv * switch_c
        dfee = - comm1 * (switch_c * delrinv + d_switch_dx*dxdr) * delr2inv

        ! Add the full non-switched coulomb pot to the softcore energy array
        sc_ener(9) = sc_ener(9) + comm1*delrinv

        ! Scaled up by oneweight
        dfee = dfee + oneweight * comm1 * delrinv * delr2inv
        if (ifcr .ne. 0) then
          b0 = (oneweight - switch_c) * delrinv
          call cr_add_dcdr_factor(i, b0*cgj)
          call cr_add_dcdr_factor(j, b0*cgi)
        end if
      else 

        ! Use electrostatic softcore potential
        delr = 1.0d0 / delrinv
        delr_n = delr**sceeorder
        if (isProcessV1) then

          ! Atoms are appearing
          denom = 1.0d0 / (delr_n + scbeta*(1.0d0-clambda))**(1.0d0/sceeorder) 
          denom_n = denom**(1+sceeorder)
          b0 = switch * denom

          ! Soft core potential goes into main electrostatic energy
          eelt = eelt + b0 * comm1

          sc_dvdl_ee = sc_dvdl_ee + &
                       switch / sceeorder * comm1 * denom_n * scbeta
          dfee = -(d_switch_dx*dxdr * comm1 * denom * delrinv) + &
                  (switch * comm1 * denom_n * delr_n * delr2inv)
          if (ifcr .ne. 0) then
            call cr_add_dcdr_factor(i, b0*cgj)
            call cr_add_dcdr_factor(j, b0*cgi)
          end if

          ! Collect lambda-dependent contributions to BAR FEP energy 
          ! differences
          if (ifmbar .ne. 0 .and. do_mbar) then
            do bar_i = 1, bar_states

              ! Remove the current lambda cont
              bar_cont(bar_i) = bar_cont(bar_i) - switch * comm1 * denom
            end do
            do bar_i = 1, bar_states
              denom = 1.0d0 / (delr_n + scbeta * &
                               (1.0d0 - bar_lambda(bar_i)))**(1.0d0/sceeorder)
              bar_cont(bar_i) = bar_cont(bar_i) + switch * comm1 * denom
            end do
          end if
        else

          ! Atoms are disappearing
          denom = 1.0d0 / (delr_n + scbeta * clambda)**(1.0d0/sceeorder)
          denom_n = denom**(1+sceeorder)
          b0 = switch * denom

          ! Soft core potential goes into main electrostatic energy
          eelt = eelt + b0 * comm1

          sc_dvdl_ee = sc_dvdl_ee - &
                       switch / sceeorder * comm1 * denom_n * scbeta
          dfee = -(d_switch_dx*dxdr * comm1 * denom * delrinv) + &
                  (switch * comm1 * denom_n * delr_n * delr2inv)
          if (ifcr .ne. 0) then
            call cr_add_dcdr_factor(i, b0*cgj)
            call cr_add_dcdr_factor(j, b0*cgi)
          end if

          ! Collect lambda-dependent contributions to
          ! Bennet Acceptance Ratio Free Energy Perturbation differences
          if (ifmbar .ne. 0 .and. do_mbar) then
            do bar_i = 1, bar_states

              ! Remove the current lambda cont
              bar_cont(bar_i) = bar_cont(bar_i) - switch * comm1 * denom
            end do
            do bar_i = 1, bar_states
              denom = 1.0d0 / (delr_n + &
                               scbeta * bar_lambda(bar_i))**(1.0d0/sceeorder)
              bar_cont(bar_i) = bar_cont(bar_i) + switch * comm1 * denom
            end do
          end if
        end if
      end if
         
      ! Inefficient, change this to use another cache later. 
      cache_r2(im_new) = cache_x(im_new) * cache_x(im_new) + &
                         cache_y(im_new) * cache_y(im_new) + &
                         cache_z(im_new) * cache_z(im_new)

      ! The following ew_directe3.h or directe4.h needs r^2 in cache_r2
      cache_df(im_new) = dfee

    end do
    ! End epilogue loop

    ! V1 uses ew_directe3.h, in which softcore atoms are treated as appearing,
    ! i.e. fully interacting at lambda=1 and 'soft' at small lambda.
    ! V0 uses ew_directe4.h, in which softcore atoms are treated as vanishing, 
    ! i.e. fully interacting at lambda=0 and 'soft' at large lambda.

    if ( isProcessV1 ) then
#     include "ew_directe3.h"
    else
#     include "ew_directe4.h"
    end if
#endif /* MPI for SOFT CORE  */
!============================================================================

#include "eedmeth2-6.h"

  end if
  ! End switch over six eedmeth options

  ! Contribute forces on the ith particle, accumulated as forces
  ! from the ith particle were added onto all other particles.
  force(1,i) = force(1,i) - dumx
  force(2,i) = force(2,i) - dumy
  force(3,i) = force(3,i) - dumz

  return

end subroutine short_ene 

#else   /* below is a drastically cut-down version, with openMP support */

!------------------------------------------------------------------------------
! get_nb_energy: the main routine for vdw, hbond, and direct space Ewald sum
!                computations. 
!------------------------------------------------------------------------------
subroutine get_nb_energy(iac, ico, ntypes, charge, cn1, cn2, cn6, force, &
                         numatoms, ipairs, ewaldcof, eedtbdns, eed_cub, &
                         eed_lin, maxnblst, eelt, evdw, ehb, dir_vir, eedvir, &
                         filter_cut, ee_type, eedmeth, dxdr, cn3, &
                         cn4, cn5, epol, dipole, field)

  use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
                     myindexlo, myindexhi, numimg, tranvec, nucgrd
  use stack
  use constants, only: zero, one, two, half, third, TWOPI, six, twelve
  use file_io_dat
  use nbips, only: teips, tvips, nnbips, rips2, ripsr, rips2r, rips6r, &
                   rips12r, aipse, aipsvc, aipsva, bipse, bipsvc, bipsva, &
                   pipsec, pipsvcc, pipsvac
!$ use omp_lib
!$ use constants, only:  omp_num_threads

  implicit none
  character(kind=1, len=13) :: routine="get_nb_energy"
#include "flocntrl.h"
#include "def_time.h"
   
  integer numatoms, maxnblst
  integer iac(*), ico(*), ntypes, ee_type, eedmeth
  _REAL_ charge(*), cn1(*), cn2(*), cn6(*)
  _REAL_ ewaldcof, eedtbdns, dxdr, eed_cub(*), eed_lin(2,*), dir_vir(3,3)
  integer ipairs(maxnblst)
  _REAL_ force(3,numatoms), eelt, epol, evdw, ehb
  _REAL_ eedvir, filter_cut, dipole(3,*), field(3,*)
  _REAL_ cn3(*), cn4(*), cn5(*)
   
  integer index, numpack, i, k, ncell_lo, ncell_hi, ntot, nvdw, nhbnd
  integer ic,j,m,n,ind,iaci,inddel
  _REAL_ del, delrinv, delr12inv
  _REAL_ b0
  _REAL_ filter_cut2
  _REAL_ comm1
  _REAL_ xktran(3,18)
  _REAL_ e3dx, e4dx, eeltl, evdwl
  _REAL_ forcel(3,numatoms)
  integer, parameter :: mask27 = 2**27 - 1
  _REAL_ delx(3), delr, delr2, cgi, delr2inv, r6, f6, f12, df, &
         dfee, dx, x, dfx(3)

!$  integer max_threads, ier

#include "../include/md.h"
  _REAL_ uips, uips2, uips4, uips2r, uips6r, uips12r
  _REAL_ pipse, dpipse, pvc, dvcu, pva, dvau
  integer itran

  _REAL_ time0, time1
   
  if (do_dir == 0) then
    return
  end if
  eelt = zero
  epol = zero
  evdw = zero
  ehb = zero
  eedvir = zero
  dir_vir(1:3,1:3) = zero
  del = one / eedtbdns
  filter_cut2 = filter_cut * filter_cut
  numpack = 1

!$ max_threads = min( 8, omp_num_threads )

  call timer_start(TIME_SHORT_ENE)

!$omp parallel private( inddel, myindexlo, myindexhi, index, &
!$omp&  k,i,ntot,nvdw,nhbnd,cgi,iaci,m,xktran, &
!$omp&  n,itran,j,delx,delr2,delrinv,x,  &
!$omp&  ind,dx,e3dx,e4dx,b0,comm1,  &
!$omp&  dfee,delr2inv,ic,r6,delr12inv,f6,f12,dfx,evdwl,eeltl,forcel)  &
!$omp&  firstprivate(numpack) num_threads(max_threads)

!$  inddel = (nucgrd-1) / OMP_GET_NUM_THREADS() + 1
!$  if (inddel == 0) inddel = 1
!$  myindexlo = 1 + OMP_GET_THREAD_NUM()*inddel
!$  myindexhi = myindexlo + inddel - 1
!$  if (OMP_GET_THREAD_NUM() == OMP_GET_NUM_THREADS()-1) myindexhi = nucgrd

  ! local versions of reductions: will be reduced in a single critical
  ! section after the end of the do loop
  evdwl = zero
  eeltl = zero
  forcel(:,:) = zero

  do index = 1,nucgrd
    if (numimg(index) <= 0) cycle
    do k = nlogrid(index), nhigrid(index)
      i = bckptr(k)
      ntot = numvdw(i) + numhbnd(i)
      if (ntot <= 0) cycle
!$    if( index < myindexlo .or. index > myindexhi ) then
!$       ! some other thread will do the work here; just need to update numpack
!$       numpack = numpack + ntot
!$       cycle
!$    endif
      nvdw = numvdw(i)
      nhbnd = numhbnd(i)

      !==========  short ene routine in-lined below:==========================
      cgi = charge(i)
      iaci = ntypes * (iac(i) - 1)
   
      xktran(1,:) = tranvec(1,:) - imagcrds(1,k)
      xktran(2,:) = tranvec(2,:) - imagcrds(2,k)
      xktran(3,:) = tranvec(3,:) - imagcrds(3,k)

      !-----------------------------------------------
      ! Loop over the 12-6 LJ and electrostatic terms 
      !-----------------------------------------------
      do m = 1, nvdw+nhbnd

        n=ipairs(m+numpack-1)
        itran=ishft(n,-27)
        n = iand(n,mask27)
        j = bckptr(n)
        delx(:) = imagcrds(:,n) + xktran(:,itran)
        delr2 = dot_product( delx, delx )

        if ( delr2 >= filter_cut2 ) cycle

        delrinv = 1.d0/sqrt(delr2)
        x = dxdr*delrinv*delr2

        ! Cubic spline on switch:
        ind = int(eedtbdns*x)
        dx = x - ind*del
        ind = 4*ind
        e3dx = dx * eed_cub(3+ind)
        e4dx = dx * dx * eed_cub(4+ind)
        b0 = delrinv*(eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
                 (e3dx + e4dx*third)*half))
        comm1 = cgi*charge(j)
        dfee = comm1*(b0 - dxdr*(eed_cub(2+ind) + e3dx + e4dx*half))
        eeltl = eeltl + comm1*b0
        delr2inv = delrinv*delrinv
        dfee = dfee*delr2inv

        ! epilogue: 12-6 LF terms

        if (m<=nvdw) then
           ic = ico(iaci+iac(j))
           r6 = delr2inv*delr2inv*delr2inv
           delr12inv = r6 * r6
           f6 = cn2(ic)*r6
           f12 = cn1(ic)*delr12inv
           evdwl = evdwl + f12 - f6
           dfee = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
        endif

        dfx(:) = delx(:)*dfee
        forcel(:,j) = forcel(:,j) + dfx(:)
        forcel(:,i) = forcel(:,i) - dfx(:)
      end do 

      !  end of what is now a single loop over j that are in the list as
      !    being close to i

      ! Contribute forces on the ith particle, accumulated as forces
      ! from the ith particle were added onto all other particles.

      !==========  short ene routine in-lined above===========================
                        
      numpack = numpack + ntot
    end do  !  k = ncell_lo,ncell_hi
  end do
!$omp critical
  force(:,:) = force(:,:) + forcel(:,:)
  evdw = evdw + evdwl
  eelt = eelt + eeltl
!$omp end critical
!$omp end parallel

  call timer_stop(TIME_SHORT_ENE)

  return

end subroutine get_nb_energy 

#endif  /* for the LES/MPI vs openMP versions */
