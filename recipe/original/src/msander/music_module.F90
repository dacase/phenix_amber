#include "../include/dprec.fh"
#include "../include/assert.fh"

! Module music handles all functionality for the
! GAL17 and vsGAL17 Pt/H2O force field
!
! Andreas Goetz
! March 7, 2017
! November 22, 2017 (MPI parallelization)
! December 1, 2017 (GAL17 parameters)
! February 2-6, 2017 (Spohr89 potential)
!
! TODO:
! - what about the virial? (serial and parallel)
! 
module music_module

  implicit none

  private
  public :: read_music_nml
  public :: print_music_settings
  public :: music_force

  ! global data from input file (music name list)
  _REAL_ :: c6_pt_o, a6_pt_o
  _REAL_ :: s_ang, r_ang, a1, a2, a3, a4, b_h_surf
  _REAL_ :: eps_gauss, bxy_gauss, bz_gauss
  integer :: n_h_surf
  logical :: debug, do_spohr89, do_dispersion, do_angle, do_theta, do_propeller, do_gauss
  character(len=2) :: pt_plane, pt_atom_type, o_atom_type, h_atom_type, vs_atom_type

  ! global data
  logical, save :: init = .false.
  logical, save :: init_plane = .false.
  _REAL_ :: surface_normal(3)
  integer :: surface_index
  _REAL_ :: surface_position
  _REAL_ :: b_gauss(3)

contains

  subroutine read_music_nml()

    use constants, only : zero, half, one, two

    implicit none

    integer, parameter :: iunit = 5  ! assume mdin is on unit 5
    integer :: ier, ifind
    namelist /music/ debug, do_spohr89, pt_atom_type, o_atom_type, h_atom_type, &
                     vs_atom_type, do_dispersion, c6_pt_o, a6_pt_o, &
                     do_angle, do_theta, do_propeller, &
                     pt_plane, s_ang, r_ang, a1, a2, a3, a4, &
                     b_h_surf, n_h_surf, &
                     do_gauss, eps_gauss, bxy_gauss, bz_gauss

    ! default namelist values
    debug = .false.
    do_spohr89 = .false.
    pt_atom_type = 'Pt'
    o_atom_type = 'OW'
    h_atom_type = 'HW'
    vs_atom_type = 'Pt'  ! use e.g. VS for Gaussian on virtual sites
    do_dispersion = .false.
    c6_pt_o = zero
    a6_pt_o = zero
    do_angle = .true.
    do_theta = .true.
    do_propeller = .true.
    pt_plane = 'xy'
    s_ang = 11.135d0     ! GAL17
    r_ang = 2.441d0      ! GAL17, Angstrom
    a1 = 15.768d0        ! GAL17, kcal/mol
    a2 = 1.594d0         ! GAL17, kcal/mol
    a3 = 1.922d0         ! GAL17, kcal/mol
    a4 = 2.838d0         ! GAL17, kcal/mol
    b_h_surf = 304.081d0 ! GAL17, kcal/mol/A^5, denoted A_{Hsurf} in paper
    n_h_surf = 5
    do_gauss = .true.
    eps_gauss = -8.901d0 ! GAL17, kcal/mol, attractive Gaussian on Pt (not VS)
    bxy_gauss = 9.331d0  ! GAL17, 1/A^2
    bz_gauss = 0.102d0   ! GAL17, 1/A^2

    ! read namelist
    rewind iunit
    call nmlsrc('music', iunit, ifind)
    if (ifind > 0) then
       rewind iunit
       read(iunit, nml=music, iostat=ier)
       REQUIRE(ier == 0)
       init = .true.
    end if

    if (do_spohr89) then
       do_dispersion = .false.
       do_angle = .false.
       do_theta = .false.
       do_propeller = .false.
       do_gauss = .false.
    end if

  end subroutine read_music_nml

  subroutine print_music_settings()

    implicit none

    if (init) then

       write(6,'(/,a)') 'MuSiC:'
       write(6,'(a)') '  (note: assumes water above platinum slab)'
       write(6,'(5x,a)')         'music namelist present in input file'
       write(6,'(5x,a,l1)')      'debug         = ', debug
       write(6,'(5x,a,l1)')      'do_spohr89       = ', do_spohr89
       write(6,'(5x,a,a)')       'pt_atom_type  = ', pt_atom_type
       write(6,'(5x,a,a)')       'vs_atom_type  = ', vs_atom_type
       write(6,'(5x,a,a)')       'o_atom_type   = ', o_atom_type
       write(6,'(5x,a,a)')       'h_atom_type   = ', h_atom_type
       write(6,'(5x,a,l1)')      'do_dispersion = ', do_dispersion
       write(6,'(5x,a,es18.11)') 'c6_pt_o       = ', c6_pt_o
       write(6,'(5x,a,es18.11)') 'a6_pt_o       = ', a6_pt_o
       write(6,'(5x,a,l1)')      'do_angle      = ', do_angle
       write(6,'(5x,a,l1)')      'do_theta      = ', do_theta
       write(6,'(5x,a,l1)')      'do_propeller  = ', do_propeller
       write(6,'(5x,a,a)')       'pt_plane      = ', pt_plane
       write(6,'(5x,a,es18.11)') 's_ang         = ', s_ang
       write(6,'(5x,a,es18.11)') 'r_ang         = ', r_ang
       write(6,'(5x,a,es18.11)') 'a1            = ', a1
       write(6,'(5x,a,es18.11)') 'a2            = ', a2
       write(6,'(5x,a,es18.11)') 'a3            = ', a3
       write(6,'(5x,a,es18.11)') 'a4            = ', a4
       write(6,'(5x,a,es18.11)') 'b_h_surf      = ', b_h_surf
       write(6,'(5x,a,i0)')      'n_h_surf      = ', n_h_surf
       write(6,'(5x,a,l1)')      'do_gauss      = ', do_gauss
       write(6,'(5x,a,es18.11)') 'eps_gauss     = ', eps_gauss
       write(6,'(5x,a,es18.11)') 'bxy_gauss     = ', bxy_gauss
       write(6,'(5x,a,es18.11)') 'bz_gauss      = ', bz_gauss

    end if
    
  end subroutine print_music_settings

#ifdef MPI
  subroutine broadcast_music_settings()

    implicit none
#include "parallel.h"
    include 'mpif.h'

    integer :: ier

    call mpi_bcast(init         , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(do_spohr89   , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(debug        , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(pt_atom_type , 2, mpi_character       , 0, commsander, ier)
    call mpi_bcast(o_atom_type  , 2, mpi_character       , 0, commsander, ier)
    call mpi_bcast(h_atom_type  , 2, mpi_character       , 0, commsander, ier)
    call mpi_bcast(vs_atom_type , 2, mpi_character       , 0, commsander, ier)
    call mpi_bcast(do_dispersion, 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(c6_pt_o      , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(a6_pt_o      , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(do_angle     , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(do_theta     , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(do_propeller , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(pt_plane     , 2, mpi_character       , 0, commsander, ier)
    call mpi_bcast(s_ang        , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(r_ang        , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(a1           , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(a2           , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(a3           , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(a4           , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(b_h_surf     , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(n_h_surf     , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(do_gauss     , 1, mpi_logical         , 0, commsander, ier)
    call mpi_bcast(eps_gauss    , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(bxy_gauss    , 1, mpi_double_precision, 0, commsander, ier)
    call mpi_bcast(bz_gauss     , 1, mpi_double_precision, 0, commsander, ier)

  end subroutine broadcast_music_settings
#endif

  subroutine music_force(ipairs, vdisp, vang, vgauss, vspohr89)

    use nblist, only : maxnblst, bckptr, imagcrds, nlogrid, nhigrid, &
                       numvdw, numhbnd, myindexlo, myindexhi, numimg
    use constants, only : zero

    implicit none 

#include "extra.h"
    integer, intent(in) :: ipairs(maxnblst)
    _REAL_, intent(out) :: vdisp, vang, vgauss, vspohr89

    integer :: numpack, ndx, ncell_lo, ncell_hi, k, i, ntot, nvdw
    _REAL_ :: xk(3)

#ifdef MPI
#include "parallel.h"
    include 'mpif.h'
    integer :: ier
    _REAL_ :: tmp_mpi
    logical, save :: first = .true.

    if (first) then
       call broadcast_music_settings()
       first = .false.
    end if
#endif

    vdisp = zero
    vang = zero
    vgauss = zero
    vspohr89 = zero
    
    ! quit if music namelist is not present
    if (.not. init) return

    if (.not. init_plane) then
       call find_plane()
       init_plane = .true.
    end if

    if (do_dispersion .or. do_gauss .or. do_spohr89) then

       ! following get_nb_energy(), see also modvdw_init()
       numpack = 1
       do ndx = myindexlo, myindexhi
          if ( numimg(ndx) > 0 ) then
             ncell_lo = nlogrid(ndx)
             ncell_hi = nhigrid(ndx)
             do k = ncell_lo, ncell_hi
                i = bckptr(k)
                xk(1:3) = imagcrds(1:3,k)
                ntot = numvdw(i) + numhbnd(i)
                nvdw = numvdw(i)
                if ( ntot > 0 ) then
                   if (do_spohr89) then
                      call spohr89_force(i, xk, ipairs(numpack), nvdw, vspohr89)
                   else
                      call dispersion_gauss_force(i, xk, ipairs(numpack), nvdw, vdisp, vgauss)
                   end if
                   numpack = numpack + ntot
                endif
             enddo
          endif
       enddo

#ifdef MPI
       call mpi_reduce(vdisp,tmp_mpi,1,mpi_double_precision,mpi_sum,0, &
            commsander,ier)
       vdisp = tmp_mpi
       call mpi_reduce(vgauss,tmp_mpi,1,mpi_double_precision,mpi_sum,0, &
            commsander,ier)
       vgauss = tmp_mpi
       call mpi_reduce(vspohr89,tmp_mpi,1,mpi_double_precision,mpi_sum,0, &
            commsander,ier)
       vspohr89 = tmp_mpi
#endif
       
       if (debug .and. master) then
          if (do_dispersion) then
             write(6,'(/,a,es18.11)') 'MuSiC Vdisp = ', vdisp
          end if
          if (do_gauss) then
             write(6,'(/,a,es18.11)') 'MuSiC Vgauss = ', vgauss
          end if
          if (do_spohr89) then
             write(6,'(/,a,es18.11)') 'MuSiC Vspohr89 = ', vspohr89
          end if
       end if

    end if

    if (do_angle) then

       call angle_force(vang)

#ifdef MPI
       call mpi_reduce(vang,tmp_mpi,1,mpi_double_precision,mpi_sum,0, &
            commsander,ier)
       vang = tmp_mpi
#endif

       if (debug .and. master) then
          write(6,'(/,a,es18.11)') 'MuSiC Vang = ', vang
       end if
       
    end if

  end subroutine music_force

  subroutine angle_force(vang)

    use memory_module, only : natom, coordinate, amber_atom_type, frc
    use constants, only : zero, half, one, two, three, four, five

    implicit none

    _REAL_, intent(inout) :: vang

    integer :: i, n, istart, iend
    _REAL_, dimension(:), pointer :: coord_o, coord_h1, coord_h2
    _REAL_ :: r_o_surf, r_h1_surf, r_h2_surf
    _REAL_ :: tmp, fermi, fermi_exp, one_fermi, fermi_grad
    _REAL_ :: dipole(3), dipnorm, dot, theta, vtheta, one_dipnorm
    _REAL_ :: vtheta_grad, acos_grad
    _REAL_ :: rn, vprop, vprop_grad

#ifdef MPI
#include "parallel.h"
    istart = iparpt(sanderrank) + 1
    iend   = iparpt(sanderrank+1) 
#else
    istart = 1
    iend   = natom
#endif
    do i = istart, iend
       
       if ( amber_atom_type(i) == o_atom_type ) then

          coord_o => coordinate(1:3,i)
          coord_h1 => coordinate(1:3,i+1)
          coord_h2 => coordinate(1:3,i+2)

          ! Fermi damping function
          r_o_surf = coord_o(surface_index) - surface_position
          tmp = r_o_surf / r_ang - one
          fermi_exp = exp(-s_ang*tmp)
          one_fermi = one / (one + fermi_exp)
          fermi = one - one_fermi

          ! gradient of Fermi 
          fermi_grad = - one_fermi * one_fermi * fermi_exp * s_ang / r_ang

          ! compute theta (cart wheel angle)
          ! this is actually the negative dipole moment
          vtheta = zero
          if (do_theta) then

             dipole = half * (coord_h1 + coord_h2) - coord_o
             dipnorm = sqrt(dot_product(dipole, dipole))
             one_dipnorm = one / dipnorm
             dot = dot_product(surface_normal, dipole)
             theta = acos(dot * one_dipnorm) ! surface normal is normalized to one
             vtheta = a1*cos(theta) + a2*cos(two*theta) &
                  + a3*cos(three*theta) + a4*cos(four*theta)

             ! force
             vtheta_grad = -a1*sin(theta) -two*a2*sin(two*theta) &
                  -three*a3*sin(three*theta) -four*a4*sin(four*theta)
             tmp = dot * one_dipnorm
             acos_grad = -one / sqrt(one - tmp*tmp)
             tmp = fermi * vtheta_grad * acos_grad * one_dipnorm
             frc(surface_index,i) = frc(surface_index,i) + tmp
             tmp = half * tmp
             frc(surface_index,i+1) = frc(surface_index,i+1) - tmp
             frc(surface_index,i+2) = frc(surface_index,i+2) - tmp
             tmp = - fermi * vtheta_grad * acos_grad * dot &
                  * one_dipnorm * one_dipnorm * one_dipnorm
             frc(1:3,i) = frc(1:3,i) + tmp * dipole(1:3)
             tmp = half * tmp
             frc(1:3,i+1) = frc(1:3,i+1) - tmp * dipole(1:3)
             frc(1:3,i+2) = frc(1:3,i+2) - tmp * dipole(1:3)

          end if

          ! correction for propeller (hydrogen repulsion)
          vprop = zero
          if (do_propeller) then

             r_h1_surf = coord_h1(surface_index) - surface_position
             rn = one
             do n = 1, n_h_surf
                rn = rn * r_h1_surf
             end do
             vprop = b_h_surf / rn
             rn = rn * r_h1_surf
             vprop_grad = -dble(n_h_surf) * b_h_surf / rn
             ! force H1
             frc(surface_index,i+1) = frc(surface_index,i+1) - vprop_grad

             r_h2_surf = coord_h2(surface_index) - surface_position
             rn = one
             do n = 1, n_h_surf
                rn = rn * r_h2_surf
             end do
             vprop = vprop + b_h_surf / rn
             rn = rn * r_h2_surf
             vprop_grad = -dble(n_h_surf) * b_h_surf / rn
             ! force H2
             frc(surface_index,i+2) = frc(surface_index,i+2) - vprop_grad

          end if

          vang = vang + fermi*vtheta + vprop

          ! force Fermi derivative
          frc(surface_index,i) = frc(surface_index,i) - fermi_grad*vtheta

       end if

    end do

  end subroutine angle_force

  subroutine find_plane()

    use memory_module, only : natom, coordinate, amber_atom_type
    use constants, only : zero, one

    implicit none

#include "extra.h"
    integer :: i
    _REAL_ :: tmp

    if (pt_plane == 'xy') then
       surface_index = 3
    else if (pt_plane == 'xz') then
       surface_index = 2
    else if (pt_plane == 'yz') then
       surface_index = 1
    else
      call sander_bomb("music: find_plane()","pt_plane is wrong", &
           "Please check your input and use xy, xz or yz.")
    end if

    surface_normal = zero
    surface_normal(surface_index) = one

    ! We assume that the surface is at the origin
    ! i.e. all atoms are at larger x, y or z values
    surface_position = -1.0d20
    do i = 1, natom
       if (amber_atom_type(i) == pt_atom_type) then
          tmp = coordinate(surface_index,i)
          if (tmp > surface_position) then
             surface_position = tmp
          end if
       end if
    end do

    if (master) then
       write(6,'(/,a,es18.11)') 'MuSiC surface_position = ', surface_position
    end if

    b_gauss(:) = bxy_gauss
    b_gauss(surface_index) = bz_gauss

  end subroutine find_plane

  subroutine dispersion_gauss_force(i, xk, ipairs, nvdw, vdisp, vgauss)

    use constants, only : one, two, six
    use nblist, only : nbfilter, tranvec, bckptr, imagcrds
    use memory_module, only : amber_atom_type, frc
    
    implicit none

    integer, intent(in) :: i, ipairs(*), nvdw
    _REAL_, intent(in) :: xk(3)
    _REAL_, intent(inout) :: vdisp, vgauss

    _REAL_ :: filter_cut2
    _REAL_ :: xktran(3,18)
    _REAL_ :: delx(3), delx2(3), r2, r4, r6, dispinv
    _REAL_ :: tmp, tmpvec(3), force(3)
    integer :: m, n, j, itran
    integer, parameter :: mask27 = 2**27 - 1
    logical :: is_pt, is_o, is_vs

    filter_cut2 = nbfilter * nbfilter

    do m = 1, 18
       xktran(1:3,m) = tranvec(1:3,m) - xk(1:3)
    end do

    do m = 1, nvdw

       n = ipairs(m)
       itran = ishft(n, -27)
       n = iand(n, mask27)
       j = bckptr(n)

       is_pt = ( amber_atom_type(j) == pt_atom_type .or. &
                 amber_atom_type(i) == pt_atom_type )
       is_o = ( amber_atom_type(j) == o_atom_type .or. &
                amber_atom_type(i) == o_atom_type )
       is_vs = ( amber_atom_type(j) == vs_atom_type .or. &
                amber_atom_type(i) == vs_atom_type )

       ! dispersion
       if ( do_dispersion .and. is_pt .and. is_o ) then

          delx(1:3) = imagcrds(1:3,n) + xktran(1:3,itran)
          r2 = delx(1)*delx(1) + delx(2)*delx(2) + delx(3)*delx(3)

          if ( r2 < filter_cut2 ) then

             r4 = r2 * r2
             r6 = r4 * r2
             dispinv = one / (r6 + a6_pt_o)
             vdisp = vdisp - c6_pt_o * dispinv

             ! force
             ! d/drk = c6_pt_o * (r6 + a6_pt_o)^(-2) * dr6/drk
             ! with dr6/drk = 6*r5*dr/drk
             ! with dr/drk = rk/r
             ! ie d/drk = 6 * c6_ptPo * r4 * delx(k) * (r6 + a6_pt_o)^(-2)
             tmp = six * c6_pt_o * r4 * dispinv * dispinv
             force(1:3) = tmp * delx(1:3)
             frc(1:3,i) = frc(1:3,i) + force(1:3)
             frc(1:3,j) = frc(1:3,j) - force(1:3)

          end if

       end if

       ! Virtual site repulsion
       if ( do_gauss .and. is_vs .and. is_o ) then

          delx(1:3) = imagcrds(1:3,n) + xktran(1:3,itran)
          delx2(1:3) = delx(1:3)*delx(1:3)
          r2 = delx2(1) + delx2(2) + delx2(3)

          if ( r2 < filter_cut2 ) then

             tmpvec(1:3) = exp(-b_gauss(1:3)*delx2(1:3))
             tmp = tmpvec(1) * tmpvec(2) * tmpvec(3)
             vgauss = vgauss + eps_gauss * tmp

             !force
             force(1:3) = -two * b_gauss(1:3) * delx(1:3) * tmp * eps_gauss
             frc(1:3,i) = frc(1:3,i) + force(1:3)
             frc(1:3,j) = frc(1:3,j) - force(1:3)

          end if

       end if

    end do

  end subroutine dispersion_gauss_force

  subroutine spohr89_force(i, xk, ipairs, nvdw, vspohr89)

    ! Implements eqs (1) through (4) of
    ! E. Spohr, J. Phys. Chem. 93 (1989) 6171-6180.
    
    use constants, only : zero, one, two
    use nblist, only : nbfilter, tranvec, bckptr, imagcrds
    use memory_module, only : amber_atom_type, frc
    
    implicit none

    integer, intent(in) :: i, ipairs(*), nvdw
    _REAL_, intent(in) :: xk(3)
    _REAL_, intent(inout) :: vspohr89

    _REAL_ :: filter_cut2
    _REAL_ :: xktran(3,18)
    _REAL_ :: delx(3), delx2(3), r, r2, rinv, rho, f_rho
    _REAL_ :: fac, tmp, tmp1, tmp2
    _REAL_ :: force(3), rho_grad(3), f_rho_grad(3)
    integer :: m, n, j, itran, ind
    integer, parameter :: mask27 = 2**27 - 1
    logical :: is_pt, is_o, is_h

    ! Spohr89 energy units are in 1.0d-19 Joule
    ! We need to convert to kcal/mol
    ! 1.0d-19 J = 1.43932643d+01 kcal/mol
    _REAL_, parameter :: to_kcal = 1.43932643d+01
    _REAL_, parameter :: small = 1.d-16 ! prevent division by zero
    _REAL_, parameter :: ah =  1.7142d0 * to_kcal
    _REAL_, parameter :: bh = -1.2777d0
    _REAL_, parameter :: ao1 = 1894.2d0 * to_kcal
    _REAL_, parameter :: ao2 = 1886.3d0 * to_kcal
    _REAL_, parameter :: ao3 = 1.0d06   * to_kcal
    _REAL_, parameter :: bo1 = -1.1004d0
    _REAL_, parameter :: bo2 = -1.0966d0
    _REAL_, parameter :: bo3 = -5.3568d0
    _REAL_, parameter :: rhofac = -0.5208d0

    filter_cut2 = nbfilter * nbfilter

    do m = 1, 18
       xktran(1:3,m) = tranvec(1:3,m) - xk(1:3)
    end do

    do m = 1, nvdw

       n = ipairs(m)
       itran = ishft(n, -27)
       n = iand(n, mask27)
       j = bckptr(n)

       is_pt = ( amber_atom_type(j) == pt_atom_type .or. &
                 amber_atom_type(i) == pt_atom_type )
       is_o = ( amber_atom_type(j) == o_atom_type .or. &
                amber_atom_type(i) == o_atom_type )
       is_h = ( amber_atom_type(j) == h_atom_type .or. &
                amber_atom_type(i) == h_atom_type )

       if ( is_pt .and. ( is_o .or. is_h ) ) then

          delx(1:3) = imagcrds(1:3,n) + xktran(1:3,itran)
          delx2(1:3) = delx(1:3)*delx(1:3)
          r2 = delx2(1) + delx2(2) + delx2(3)

          if ( r2 < filter_cut2 ) then

             r = sqrt(r2)
             rinv = one / r

             if ( is_o ) then
             
                ! rho: projection of distance vector onto surface plane
                rho = zero
                do ind = 1, 3
                   fac = one
                   if ( ind == surface_index) then
                      fac = zero
                   end if
                   rho = rho + fac*delx2(ind)
                end do
                rho = sqrt(rho)
                if (rho < small) then    ! in this case delx is also zero
                   rho_grad(1:3) = zero  ! prevent division by zero
                else
                   rho_grad (1:3) = delx(1:3) / rho
                end if
                rho_grad(surface_index) = zero

                f_rho = exp(rhofac*rho*rho)
                tmp = f_rho * rhofac * two * rho
                f_rho_grad(1:3) = tmp * rho_grad(1:3)

                tmp1 = ao1*exp(bo1*r)
                tmp2 = ao2*exp(bo2*r)
                tmp = f_rho * (tmp1 - tmp2)
                force(1:3) = f_rho * (tmp1*bo1 - tmp2*bo2) * rinv * delx(1:3)
                force(1:3) = force(1:3) + (tmp1 - tmp2) * f_rho_grad(1:3)
                
                tmp1 = ao3*exp(bo3*r)
                tmp = tmp + tmp1*(one-f_rho)
                force(1:3) = force(1:3) + (one-f_rho)*tmp1*bo3 * rinv * delx(1:3)
                force(1:3) = force(1:3) - tmp1 * f_rho_grad(1:3)

                frc(1:3,i) = frc(1:3,i) + force(1:3)
                frc(1:3,j) = frc(1:3,j) - force(1:3)

             else  ! must be h

                tmp = ah*exp(bh*r)

                !force
                force(1:3) = tmp * bh * rinv * delx(1:3)
                frc(1:3,i) = frc(1:3,i) + force(1:3)
                frc(1:3,j) = frc(1:3,j) - force(1:3)

             end if

             vspohr89 = vspohr89 + tmp

          end if

       end if

    end do

  end subroutine spohr89_force
  
end module music_module
