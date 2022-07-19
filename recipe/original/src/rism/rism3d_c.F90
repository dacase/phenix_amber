!<compile=optimized>
#include "../include/dprec.fh"

!> 3D-RISM solver.
!! This defines the 3D-RISM type and associated subroutines.  All type
!! elements are public.  In general, read but do not write these
!! variables.  This provides an object-orientiented interface without
!! needing a function to access every variable.
!!
!! Features of this solver include:
!! o Multiple closures w/ temperature derivatives: KH, HNC, PSE-n
!! o Temperature derivative expressed as T*d/dT
!! o MDIIS accelerated solutions
!! o Optional cutoffs
!! o Analytic forces
!! o Variable grid size and dynamic memory allocation
!! o MPI support
!! o units:  energy       [kT]
!!           distances    [A]         (Angstroms)
!!           site charges [sqrt(kT A)]
!!           temperature  [K]
!!           density      [#/A^3]
!!           mass         [au]
!! o To convert [e] to [sqrt(kT A)] * sqrt(COULOMB_CONST_E/ KB / temperature)

module rism3d_c
  use rism3d_solute_c
  use rism3d_solvent_c
  use rism3d_potential_c
  use rism3d_grid_c
  use rism3d_closure_c
  use rism_report_c
  use mdiis_c
  use rism3d_fft_c

  use rism3d_opendx
#ifdef RISM3D_DEBUG
  !    use rism3d_debug_c
#endif
  implicit none
#include "def_time.h"

  type rism3d
     !! Solute/solvent information.

     !> Solute object.
     type(rism3d_solute) :: solute
     !> Solvent object.
     type(rism3d_solvent) :: solvent
     !> Potential object.
     type(rism3d_potential) :: potential
     !> Grid object.
     type(rism3d_grid) :: grid
     !> Closure object.
     type(rism3d_closure) :: closure

     !> List of closure names to use in order.  Only the last closure
     !! is used for thermodynamic output.  This can be used to
     !! progressively increase the order of the closure to aid
     !! convergence.
     character(len = 8), pointer :: closureList(:) => NULL()

     ! TIMERS.  Subtimers only account for computation.  We ignore setup etc.
     ! timer :: timer for this class.  Activated for all public routines
     ! resizeTimer :: time to resize solvent box
     ! reorientTimer :: time to reorient solute
     ! cuvpropTimer :: time to propagate Cuv solution
     ! fftTimer :: specifically times FFT calculation
     ! solveTimer :: specifically times rism1d_solve calculation
     ! solve3DRISMTimer :: specifically times solve3DRISM calculation
     ! single3DRISMsolutionTimer :: specifically times single3DRISMsolution calculation
     ! thermoTimer :: specifically times thermodynamics calculations
     ! forceTimer :: specifically times force calculation
     ! excessChemicalPotentialTimer :: specificall times excess chemical potential calculation

     !! private !(should be)

     ! FFTW options

     ! FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE
     integer :: fftw_planner = FFT_PATIENT
     ! .true.  - use aligned memory and to enable SIMD;
     ! .false. - don't use aligned memory
     logical :: fft_aligned = .true.
     ! Transpose site number and spatial data locally before and after FFT.
     logical :: fftw_localtrans = .true.

     !> Output verbosity.  Useful for debugging.
     !! 0 - no ouput
     !! 1 - memory allocation and steps for convergence
     !! 2 - 1 + convergence progress
     integer :: verbose = 0

     ! This is a bit ugly and there may be a better solution.  We
     ! need to keep track of the number of solutions for both charged
     ! and un-charged solutes.  When we change between the solutes we
     ! set the nsolutions pointer to the appropriate variable.
     ! However, the 'target' attribute is not allowed in type
     ! definitions so these variables have to be pointers and we have
     ! to allocate memory for them.

     !> Number of times full solutions have been calculated.
     integer, pointer :: nsolution => NULL()
     !> Number of times full solutions with a charged solute have been
     !! calculated.
     integer, pointer :: nsolutionChg => NULL()
     !> Number of times full solutions with an uncharged solute have
     !! been calculated.
     integer, pointer :: nsolutionNoChg => NULL()

     !> Number of past direct correlation function time step saves.
     integer :: ncuvsteps ! numDCFsteps

     !> Fixed box size for 3D-RISM.
     _REAL_ :: fixedBoxDimensionsR(3)
     !> Number of Cartesian grid points in each dimension for a fixed box size.
     integer :: fixedNumGridPoints(3)

     !> Number of vectors used for MDIIS (consequently, the number of
     !! copies of CUV we need to keep for MDIIS).
     integer :: NVec
     !> MDIIS implementation to use.
     integer :: mdiis_method
     type(mdiis) :: mdiis_o

     !> 'Step size' for MDIIS.
     _REAL_ :: deloz = 0.7d0
     !> Restart threshold factor. Ratio of the current residual to the
     !! minimum residual in the basis that causes a restart.
     _REAL_ :: mdiis_restart

     !! MPI Support !!
     integer :: mpirank = 0, mpicomm = 0, mpisize = 1

     !! LARGE ARRAYS !!
     !
     ! all arrays are declared as pointers to ensure we can reallocate them as necessary
     !

     ! xvva       :: solvent chi interpolated for our grid size
     ! guv        :: solvent distribution function
     ! huv        :: guv - 1
     ! cuv        :: solvent direct correlation function and points to the
     !              current active solution in cuvWRK
     ! cuvres     :: residual value for cuv calculation and points to the
     !              current active solution in cuvresWRK.
     ! cuvWRK     :: Working Cuv memory.  Holds Cuv from previous iterations.
     ! cuvresWRK  :: Working Cuvres memory.  Holds Cuvres from previous iterations.
     ! oldcuv     :: previous solutions of cuv. Points to oldcuvChg
     !              or oldcuvNoChg depending on the charge state of the
     !              calculation.
     ! oldcuvChg  :: previous solutions for the standard charged system
     ! oldcuvNoChg :: previous solutions for the chargeless system.  This is only allocated
     !                if _unsetCharges() is called
     _REAL_, pointer :: xvva(:) => NULL(), &
          oldcuv(:, :, :, :, :) => NULL(), &
          oldcuvChg(:, :, :, :, :) => NULL(), &
          oldcuvNoChg(:, :, :, :, :) => NULL(), &
          cuv(:, :, :, :) => NULL(), cuvWRK(:, :, :, :, :) => NULL(), &
          cuvres(:, :) => NULL(), cuvresWRK(:, :, :) => NULL()


     _REAL_, pointer :: guv(:, :) => NULL(), huv(:, :) => NULL()

     ! cuvk        :: k-space Cuv solution from 3D-RISM
     !               solution. NOTE: we should consider using Huv or
     !               Guv memory instead.  However, it has to be
     !               checked first that it is not used for any thermodynamics calculations

     _REAL_, pointer :: cuvk(:, :) => NULL()

     ! fft :: fft object for standard 3D-RISM solution
     type(rism3d_fft) :: fft

     !> If true, a periodic 3D-RISM calculation is performed. This
     !! primarily differs from infinite dilution 3D-RISM by using
     !! Ewald sum potential in place of Coulombic potential and
     !! invoking the minimum image convention while calculating both
     !! the Ewald sum and Lennard-Jones potentials.
     logical :: periodic = .false.

     !> Lengths and interior angles of the unit cell. For aperiodic
     !! systems, the interior angles are always 90 degrees.
     _REAL_ :: unitCellDimensions(6)
     
     !> The abbreviated label of the periodic potential function used
     !! for periodic calculations. See rism3d_potential for valid values.
     character(len=255) :: periodicPotential = ""

  end type rism3d

  public :: rism3d_new, rism3d_destroy, rism3d_calculateSolution, rism3d_force, &
       rism3d_excessChemicalPotential_tot, rism3d_excessChemicalPotential, &
       rism3d_setclosure, rism3d_setverbosity, rism3d_setcut, rism3d_setmdiis

  private :: resizeBox, reallocateBox, &
       interpolateSolventSusceptibility, &
       solve3DRISM, single3DRISMsolution,  &
       guessDCF, updateDCFguessHistory

contains


  !> Constructor - precalculates the solute solvent terms that are not
  !! configuration dependent and sets box parameters.
  !!
  !! The unit cell parameters always give the size of the box. 
  !!
  !! If grdspc(1:3) is set, this array gives an approximate grid
  !! spacing.  The actual spacing will be set that an exact number of 
  !! grids spans the box, and so that the number of grid points is even.  
  !! (In addtion, for MPI runs, the number of grids along y and z will 
  !! be adjusted to be a multiple of the number of MPI threads.)
  !!
  !! Alternatively, if ng3(1:3) is set, these values will be used for
  !! the number of grids in each direction, and grdspc() will be ignored.
  !!
  !! If this is an MPI run, supply the MPI communicator.  Only the rank
  !! 0 parameters will be used. However, due to the limitations of
  !! pre-Fortran2003, the closure must be the same length on all
  !! processes. The values and number of elements for the closure list
  !! on non-rank 0 processes still do not matter.
  !!
  !! IN:
  !!   this :: new rism3d object
  !!   solu :: 3D-RISM solute object
  !!   solv :: 3D-RISM solvent object
  !!   ncuvsteps :: number of past cuv time steps saves
  !!   closure :: list of closures. Closures may be KH, HNC or PSEn
  !!              where n is an integer. Ensure the length attribute is
  !!              the same on all processes.
  !!   cut     :: distance cutoff for potential and force calculations
  !!   mdiis_nvec :: number of MDIIS vectors (previous iterations) to keep
  !!   mdiis_del :: scaling factor applied to estimated gradient (residual)
  !!   mdiis_method :: which implementation of the algorithm
  !!   chargeSmear :: Charge smearing parameter for long-range
  !!       asymtotics and Ewald, typically eta in the literature
  !!   o_grdspc :: (optional) grid spacing for the solvent box in each dimension
  !!   o_ng3    :: (optional) number of grid points in each dimension
  !!   o_mpicomm :: (optional) MPI communicator
  !!   o_periodic :: (optional) periodic electric potential to use, if any
  !!   o_unitCellDimensions :: (optional) geometry of the system unit cell

  subroutine rism3d_new(this, solute, solvent, ncuvsteps, &
       closure, cut, mdiis_nvec, mdiis_del, mdiis_method, mdiis_restart, &
       chargeSmear, o_grdspc, o_ng3, o_mpicomm, &
       o_periodic, o_unitCellDimensions)
    use rism3d_solute_c
    use rism3d_solvent_c
    use safemem
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(rism3d), intent(inout) :: this
    type(rism3d_solute), intent(in), target :: solute
    type(rism3d_solvent), intent(in), target :: solvent
    integer, intent(in) :: ncuvsteps
    character(len = *), intent(in) :: closure(:)
    _REAL_, intent(in) :: cut
    integer, intent(in) :: mdiis_nvec, mdiis_method
    _REAL_, intent(in) :: mdiis_del, mdiis_restart
    _REAL_, optional, intent(in) :: o_grdspc(3)
    integer, optional, intent(in) :: o_ng3(3)
    integer, optional, intent(in) :: o_mpicomm
    character(len = *), optional, intent(in) :: o_periodic
    _REAL_, optional, intent(in) :: o_unitCellDimensions(6)
    _REAL_, intent(in) :: chargeSmear
    ! temporary copies
    character(len = len(closure)), pointer :: t_closure(:)
    _REAL_ :: t_cut
    integer :: t_mdiis_nvec, t_mdiis_method
    _REAL_ :: t_mdiis_del, t_mdiis_restart
    _REAL_ :: t_grdspc(3)
    integer :: t_ng3(3)
    _REAL_ :: t_unitCellDimensions(6)
    integer :: t_mpicomm
    _REAL_ :: t_chargeSmear
    integer :: nclosure
    integer :: err

    ! MPI set up starts by obtaining rank and size.  Temporary copies
    ! of input parameters that do not directly set object variables
    ! are made.  These are they broadcast to the rank > 0 processes.
    ! Then all processes complete the intitialization procedure
    ! using the temporary copies.  This leaves the input parameters
    ! untouched.

    nullify(t_closure)

    if (present(o_periodic)) then
       if (o_periodic /= '') then
          this%periodic = .true.
          this%periodicPotential = o_periodic
       end if
    end if

    ! GET RANK AND SIZE
    this%mpicomm = 0
    this%mpisize = 1
    this%mpirank = 0
#ifdef MPI
    if (present(o_mpicomm)) then
       this%mpicomm = o_mpicomm
       if (this%mpicomm == MPI_COMM_NULL) &
            call rism_report_error("RISM3D: received NULL MPI communicator")
       call mpi_comm_rank(this%mpicomm, this%mpirank, err)
       if (err /= 0) call rism_report_error &
            ("(a,i8)", "RISM3D: could not get MPI rank for communicator ", this%mpicomm)
       call mpi_comm_size(this%mpicomm, this%mpisize, err)
       if (err /= 0) call rism_report_error &
            ("(a,i8)", "RISM3D: could not get MPI size for communicator ", this%mpicomm)
       call rism_report_mpi(this%mpicomm)
    end if
#endif /*MPI*/
    ! MAKE TEMPORARY COPIES
    if (this%mpirank == 0) then

       call rism3d_solute_clone(solute, this%solute)
       call rism3d_solvent_clone(solvent, this%solvent)
       this%ncuvsteps = ncuvsteps
       nclosure = size(closure)
       t_closure => safemem_realloc(t_closure, len(closure), nclosure)
       t_closure = closure
       t_cut = cut
       t_mdiis_nvec = mdiis_nvec
       t_mdiis_method = mdiis_method
       t_mdiis_del = mdiis_del
       t_mdiis_restart = mdiis_restart
       ! check box parameters
       if (present(o_grdspc)) then
          t_grdspc = o_grdspc
          if (present(o_ng3)) &
             call rism_report_error("RISM3D: do not set both GRDSPC and NG3")
       else if (present(o_ng3)) then
          t_ng3 = o_ng3
          if (present(o_grdspc)) &
             call rism_report_error( "RISM3D: do not set both GRDSPC and NG3")
       else
          call rism_report_error( "RISM3D: must set either GRDSPC or NG3")
       end if
       if (present(o_unitCellDimensions)) then
          t_unitCellDimensions = o_unitCellDimensions
       end if
       t_chargeSmear = chargeSmear
    end if
#ifdef MPI
    ! BROADCAST PARAMETERS
    ! set solu on all processes
    call rism3d_solute_mpi_clone(this%solute, this%mpirank, this%mpicomm)
    ! set solv on all processes
    call rism3d_solvent_mpi_clone(this%solvent, this%mpirank, this%mpicomm)
    ! set ncuvstpes on all processes
    call mpi_bcast(this%ncuvsteps, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /= 0) call rism_report_error("RISM3D: broadcast NCUVSTEPS in constructor failed")
    call mpi_bcast(nclosure, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error &
         ("RISM3D interface: could not broadcast PROGRESS")
    if (this%mpirank/=0) &
         t_closure => safemem_realloc(t_closure, len(t_closure), nclosure)
    call mpi_bcast(t_closure, len(t_closure) * nclosure, mpi_character, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast CLOSURE in constructor failed")
    call mpi_bcast(t_cut, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast CUT in constructor failed")
    call mpi_bcast(t_mdiis_nvec, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_NVEC in constructor failed")
    call mpi_bcast(t_mdiis_del, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_DEL in constructor failed")
    call mpi_bcast(t_mdiis_restart, 1, mpi_double_precision, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_RESTART in constructor failed")
    call mpi_bcast(t_mdiis_method, 1, mpi_integer, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast MDIIS_METHOD in constructor failed")
    if (present(o_grdspc)) then
       call mpi_bcast(t_grdspc, 3, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast GRDSPC in constructor failed")
    else
       call mpi_bcast(t_ng3, 3, mpi_integer, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast NG3 in constructor failed")
    end if
    if (present(o_unitCellDimensions)) then
       call mpi_bcast(t_unitCellDimensions, 6, mpi_double_precision, 0, this%mpicomm, err)
       if (err /=0) call rism_report_error("RISM3D: broadcast UNITCELLDIMENSIONS in constructor failed")
    end if
    call mpi_bcast(t_chargeSmear, 1, mpi_double, 0, this%mpicomm, err)
    if (err /=0) call rism_report_error("RISM3D: broadcast chargeSmear in constructor failed")
#endif /*MPI*/
    ! INITIALIZE
    call rism3d_grid_new(this%grid, this%mpicomm)
    call rism3d_setmdiis(this, t_mdiis_nvec, t_mdiis_del, t_mdiis_method, t_mdiis_restart)
    call rism3d_potential_new(this%potential, this%grid, this%solvent, this%solute, 0d0, &
         this%fft, this%periodicPotential, chargeSmear)

#ifdef MPI
    call mdiis_new_mpi(this%mdiis_o, this%mdiis_method, &
         this%deloz, 0d0, &
         this%MDIIS_restart, &
         this%mpirank, this%mpisize, this%mpicomm)
#else
    call mdiis_new(this%mdiis_o, this%mdiis_method, &
         this%deloz, 0d0, &
         this%MDIIS_restart)
#endif /*MPI*/

    call rism3d_setcut(this, t_cut)
    call rism3d_setclosurelist(this, t_closure)
    if (present(o_grdspc)) then
       call rism3d_grid_setSpacing(this%grid, t_grdspc)
    else
       this%fixedNumGridPoints = t_ng3
    end if

    if (present(o_unitCellDimensions)) then
       this%unitCellDimensions = t_unitCellDimensions
       call rism3d_grid_setUnitCellDimensions(this%grid, this%unitCellDimensions, this%periodic)
    end if

    allocate(this%nsolutionChg, this%nsolutionNoChg)
    this%nsolutionChg = 0
    this%nsolutionNoChg = 0
    this%nsolution => this%nsolutionChg

    
    call rism3d_fft_global_init()
    this%fftw_planner = FFT_MEASURE
#if defined(MPI)
    this%fft_aligned = .false.
#else
    this%fft_aligned = .true.
#endif
    this%fftw_localtrans = .true.

    ! Clean up locally allocated temporary memory.
    if (safemem_dealloc(t_closure) /= 0) &
         call rism_report_error("RISM3D:NEW: failed to deallocate t_closure")
#ifdef RISM3D_DEBUG
    call rism3d_debug_new(this%grid, this%solvent, this%mpirank, this%mpisize, this%mpicomm)
#endif

  end subroutine rism3d_new

  !> Sets the closure list and sets the current closure to the first one
  !! in the list.  When there is no previous solution to work from, the
  !! solver will use each closure in the list in turn. By choosing the
  !! list to increase in order, it makes it possible to converge
  !! otherwise difficult closures. Only the last closure is used for
  !! thermodynamic output.
  !! IN:
  !!   this :: rism3d object
  !!   closure :: array of closure types (see closure enumeration).
  subroutine rism3d_setclosurelist(this, closure)
    implicit none
    type(rism3d), intent(inout) :: this
    character(len = *), intent(in) :: closure(:)
    this%closureList => safemem_realloc(this%closureList, len(this%closureList), &
         ubound(closure, 1))
    this%closureList = closure
    call rism3d_setclosure(this, this%closureList(1))
  end subroutine rism3d_setclosurelist


  !> Sets the closure type.
  !! IN:
  !!   this :: rism3d object
  !!   closure :: closure type (see closure enumeration).
  subroutine rism3d_setclosure(this, closure)
    implicit none
    type(rism3d), intent(inout) :: this
    character(len = *), intent(in) :: closure
    call rism3d_closure_destroy(this%closure)
    call rism3d_closure_new(this%closure, closure, this%potential)
  end subroutine rism3d_setclosure


  !> Sets verbosity of output.
  !! IN:
  !!   this :: rism3d object
  !!   verbosity :: 0 - no output
  !!                1 - memory allocation and steps for convergence
  !!                2 - 1 + convergence progress
  subroutine rism3d_setverbosity(this, verbosity)
    implicit none
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: verbosity
    this%verbose = verbosity
  end subroutine rism3d_setverbosity


  !> Sets the cut off distance for periodic potential and force calculations.
  !! IN:
  !!   this :: rism3d object
  !!   cut     :: distance cutoff for potential and force calculations
  subroutine rism3d_setcut(this, cut)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: cut
    call rism3d_potential_setCut_ljdistance(this%potential, cut)
  end subroutine rism3d_setcut


  !> Sets MDIIS parameters
  !! IN:
  !!   this :: rism3d object!
  !!   nvec :: number of MDIIS vectors (previous iterations) to keep
  !!   del :: scaling factor (step size) applied to estimated gradient (residual)
  !!   method :: which implementation of the algorithm
  !!   restart :: restart threshold factor. Ratio of the current residual to the
  !!              minimum residual in the basis that causes a restart
  subroutine rism3d_setmdiis(this, nvec, del, method, restart)
    implicit none
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: nvec, method
    _REAL_, intent(in) :: del, restart
    this%NVec = nvec
    this%deloz = del
    this%mdiis_method = method
    this%mdiis_restart = restart
  end subroutine rism3d_setmdiis


  !> Sets solute coordinates.
  !! IN:
  !!   this :: rism3d object
  !!   ratu :: coordinates
  subroutine rism3d_setCoord(this, solutePositions)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: solutePositions(:, :)
    call rism3d_solute_setCoord(this%solute, solutePositions)
  end subroutine rism3d_setCoord


  !> Sets all solute partial charges to zero, resets MDIIS and wipes out
  !! working memory.
  !! IN:
  !!   this :: rism3d object
  subroutine rism3d_unsetCharges(this)
    implicit none
    type(rism3d), intent(inout) :: this
    integer :: i
    ! reset MDIIS.  This makes the working vector index 1
    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    ! turn off charges
    call rism3d_solute_unsetCharges(this%solute)
    ! Use the number of no charge solutions
    this%nsolution => this%nsolutionNoChg
    ! Use no charge previous soluitions
    this%oldcuv => this%oldcuvNoChg
    ! Attempt to restore last solution without charge here.  If this
    ! fails, attempt again at the end of reallocateBox
    if (this%nsolutionNoChg > 0) then
       if (all(ubound(this%oldcuv(:,:,:,:,1)) .eq. ubound(this%cuv))) then
          call dcopy(product(ubound(this%cuv)), this%oldcuv, 1,this%cuv, 1)
       end if
    end if
  end subroutine rism3d_unsetCharges


  !> Sets all solute partial charges to to their original
  !! values. (Undoes rism3d_unsetCharge().)
  !! IN:
  !!   this :: rism3d object
  subroutine rism3d_resetCharges(this)
    implicit none
    type(rism3d), intent(inout) :: this
    integer :: i
    ! reset MDIIS.  This makes the working vector index 1
    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    ! get back the charges
    call rism3d_solute_resetCharges(this%solute)
    ! restore the number of previous solutions
    this%nsolution => this%nsolutionChg
    ! point to previous charged solutions
    this%oldcuv => this%oldcuvChg
    ! Attempt to restore last solution with charge here.  If this
    ! fails, attempt again at the end of reallocateBox
    if (this%nsolutionNoChg > 0) then
       if (all(ubound(this%oldcuv(:,:,:,:,1)) .eq. ubound(this%cuv))) then
          call dcopy(product(ubound(this%cuv)), this%oldcuv, 1,this%cuv, 1)
       end if
    end if
  end subroutine rism3d_resetCharges


  !> Calculates the full 3D-RISM solvent distribution.  This is required to
  !! calculate thermodynamic quantities.
  !! @param[in,out] this rism3d object.
  !! @param[in,out] ksave Save intermediate results every ksave
  !!            interations (0 means no saves).
  !! @param[in] kshow Print parameter for relaxation steps every kshow
  !!            iteration (0 means no print).
  !! @param[in] maxSteps Maximum number of rism relaxation steps.
  !! @param[in] tolerance Convergence tolerances. There should be one
  !!          tolerance per closure in the closure list.
  subroutine rism3d_calculateSolution(this, ksave, kshow, maxSteps, &
          tolerance, ng3, verbose)
    use constants_rism, only : pi
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: ksave, kshow, maxSteps
    _REAL_, intent(in) :: tolerance(:)
    integer, intent(in) :: ng3(3), verbose

    _REAL_ :: com(3)
    ! iclosure :: counter for closures
    integer :: iclosure

    _REAL_ :: offset(3)
    integer :: id, iu

    ! 1) Quick check that the tolerance list is of the correct length.
    if (ubound(tolerance, 1) /= ubound(this%closureList, 1)) &
         call rism_report_error("(a,i3,a,i3)", &
         "RISM3D_SOLVE: number of tolerances, ", &
         ubound(tolerance, 1), ", is not equal to numer of closures, ", &
         ubound(this%closureList, 1))

    ! 2) Get the box size (assumed to be fixed here)
    if (this%nsolution == 0) then
       call timer_start(TIME_RESIZE)
       call resizeBox(this,ng3)
       call timer_stop(TIME_RESIZE)

       if(verbose >= 0 ) then
         call rism_report_message("||Setting solvation box to")
         call rism_report_message("(3(a,i10))", "|grid size: ", &
            this%grid%globalDimsR(1), " X ", this%grid%globalDimsR(2), &
            " X ", this%grid%globalDimsR(3))
         call rism_report_message("(3(a,f10.3))", "|box size [A]:  ", &
            this%grid%boxLength(1), " X ", this%grid%boxLength(2), &
            " X ", this%grid%boxLength(3))
         call rism_report_message("(3(a,f10.3))", "|grid spacing [A]: ", &
            this%grid%spacing(1), " X ", this%grid%spacing(2), &
            " X ", this%grid%spacing(3))
         call rism_report_message("(3(a,f10.3))", "|internal angles [°]:  ", &
            this%grid%unitCellAngles(1) * 180 / pi, ", ", &
            this%grid%unitCellAngles(2) * 180 / pi, ", ", &
            this%grid%unitCellAngles(3) * 180 / pi)
         call rism_report_message("(a,f10.3)", "|inscribed sphere radius [A]: ",&
            this%grid%inscribedSphereRadius)
         call flush(rism_report_getmunit())
       end if
    end if

    ! 2a) Check what kind of information is in the xvv file:
    call check_xvv_info( this%potential )
    
    ! 3) Calculate electrostatic and Lennard-Jones potential about the
    ! solute.
    call rism3d_potential_calc(this%potential)

    ! 4) Propagate previously saved solute-solvent DCF solutions to
    ! create an initial guess for this solution.
    call timer_start(TIME_CUVPROP)
    call guessDCF(this)
    call timer_stop(TIME_CUVPROP)


    ! 5) Calculate 3D-RISM solution using MDIIS.
    ! If the user has to provide a list of closures, use it only if
    ! this is the first solution (nsolution == 0) or solution
    ! propagation is turned off (ncuvsteps == 0). Otherwise, the
    ! current closure will be the last one in the list.
    if (this%nsolution == 0 .or. this%ncuvsteps == 0) then
       do iclosure = 1, size(this%closureList)
          if (this%verbose >= 1) &
               call rism_report_message("|Switching to "// &
               trim(this%closureList(iclosure))//" closure")
          call rism3d_setClosure(this, this%closureList(iclosure))
          call timer_start(TIME_RXRISM)
          call solve3DRISM(this, ksave, kshow, maxSteps, &
               tolerance(iclosure))
          call timer_stop(TIME_RXRISM)
          ! Increment nsolution and ncuvsteps to ensure the previous
          ! closure solution is used.
          if (iclosure == 1) then
             this%nsolution = this%nsolution + 1
             this%ncuvsteps = this%ncuvsteps + 1
          end if
       end do
       this%nsolution = this%nsolution - 1
       this%ncuvsteps = this%ncuvsteps - 1
    else
       call timer_start(TIME_RXRISM)
       call solve3DRISM(this, ksave, kshow, maxSteps, &
            tolerance(size(tolerance)))
       call timer_stop(TIME_RXRISM)
    end if

    ! 6) Update stored variables.
    call timer_start(TIME_CUVPROP)
    this%nsolution = this%nsolution + 1
    call updateDCFguessHistory(this)
    call timer_stop(TIME_CUVPROP)

  end subroutine rism3d_calculateSolution

  !> Calculates the forces on the solute contributed by the solvent according
  !! to 3D-RISM.  Just a wrapper for rism3d_closure_force().
  !! IN:
  !!   this :: rism3d object with computed solution
  !!   ff   :: 3D-RISM forces
  subroutine rism3d_force(this, ff)
    implicit none
    type(rism3d):: this
    _REAL_, intent(out) :: ff(3, this%solute%numAtoms)
    call timer_start(TIME_FF)
    call rism3d_closure_force(this%closure, ff, this%guv, this%periodicPotential)
    call timer_stop(TIME_FF)
    !!!!!!!
    !! - uncomment this to have a proper force check
    !! call rism3d_checkForceNumDeriv(this%closure, ff, this%guv)
    !!!!!!
  end subroutine rism3d_force


  !> Calculate the excess chemical potential for each solvent species
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    excess chemical potential of solvation for each solvent species
  function rism3d_excessChemicalPotential(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)

    excessChemicalPotential = &
       rism3d_closure_excessChemicalPotential(this%closure, &
       this%huv, this%cuv(:, :, :, :))
  end function rism3d_excessChemicalPotential

  !> Calculate the total excess chemical potential of solvation
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    total excess chemical potential of solvation
  function rism3d_excessChemicalPotential_tot(this) result(excessChemicalPotential)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: excessChemicalPotential

    excessChemicalPotential = sum(rism3d_excessChemicalPotential(this))
  end function rism3d_excessChemicalPotential_tot

  !> Calculate the solvation interaction energy: de = density sum g*u for
  !! each solvent site.  I.e., the direct intection potential energy of
  !! solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    the contribution of each solvent site to the total solvation interaction
  !!    energy [kT]
  function rism3d_solventPotEne(this) result(ene)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: ene(this%solvent%numAtomTypes)
    ene = rism3d_closure_solvPotEne(this%closure, this%guv)
  end function rism3d_solventPotEne


  !> Calculate the total solvation interaction energy: de = density sum g*u for
  !! each solvent site.  I.e., the direct intection potential energy of
  !! solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    the total solvent-solute potential energy [kT]
  function rism3d_solventPotEne_tot(this) result(ene)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: ene
    ene = sum(rism3d_solventPotEne(this))
  end function rism3d_solventPotEne_tot

  !! Calculating the partial molar volume of solute.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!   partial molar volume
  function rism3d_partialMolarVolume(this) result(partialMolarVolume)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: partialMolarVolume
    partialMolarVolume = rism3d_closure_partialMolarVolume(this%closure, this%cuv(:, :, :, :))
  end function rism3d_partialMolarVolume

  !> Calculating excess number of each solvent type associated with
  !! the solute.
  !! @param[in,out] this rism3d object with computed solution.
  !! @return Excess number of each solvent type associated with the solute.
  function rism3d_excessParticles(this) result(num)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: num(this%solvent%numAtomTypes)

    num = rism3d_closure_excessParticles(this%closure, this%guv)

  end function rism3d_excessParticles

  !! Calculate the Kirkwood-Buff integral for the solute. This is the
  !! all space integral of huv.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! IN:
  !!    this :: rism3d object with computed solution
  !! OUT:
  !!    Kirkwood-Buff integeral for each solvent site
  function rism3d_kirkwoodBuff(this) result(kb)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: kb(this%solvent%numAtomTypes)

    kb = rism3d_closure_kirkwoodBuff(this%closure, this%guv)
  end function rism3d_kirkwoodBuff

  !! Calculates the direct correlation function integral for the solute. This is the
  !! all space integral of cuv.
  !! IN:
  !!   this :: rism3d object with computed solution
  !! OUT:
  !!    DCF integeral for each solvent site
  function rism3d_DCFintegral(this) result(dcfi)
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_ :: dcfi(this%solvent%numAtomTypes)
    dcfi = rism3d_closure_DCFintegral(this%closure, this%cuv)
  end function rism3d_DCFintegral

!!!!!
  !! DEALLOCATE
!!!!!

!!!!!
  !! deconstructor - frees all memory
!!!!!
  subroutine rism3d_destroy(this)
    use safemem
    implicit none
    type(rism3d) :: this

    call rism3d_solvent_destroy(this%solvent)
    call rism3d_solute_destroy(this%solute)
    call rism3d_potential_destroy(this%potential)
    call rism3d_grid_destroy(this%grid)
    call rism3d_closure_destroy(this%closure)
    call mdiis_destroy(this%mdiis_o)

    if (safemem_dealloc(this%xvva) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate XVVA")
    if (safemem_dealloc(this%cuvWRK) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVWRK")
    if (safemem_dealloc(this%oldcuvChg) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate OLDCUVCHG")
    if (safemem_dealloc(this%oldcuvNoChg) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate OLDCUVNOCHG")
    if (safemem_dealloc(this%cuvresWRK) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVRESWRK")

    if (safemem_dealloc(this%closureList) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CLOSURELIST")
    if (associated(this%nsolutionChg)) &
         deallocate(this%nsolutionChg)
    if (associated(this%nsolutionNoChg)) &
         deallocate(this%nsolutionNoChg)
    nullify(this%cuv)
    nullify(this%cuvres)
    nullify(this%oldcuv)
    nullify(this%nsolution)

    if (safemem_dealloc(this%guv, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate GUV")
    if (safemem_dealloc(this%huv, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate HUV")
    if (safemem_dealloc(this%cuvk, o_aligned = .true.) /= 0) &
         call rism_report_error("RISM3D: failed to deallocate CUVK")
    call rism3d_fft_destroy(this%fft)
    call rism3d_fft_global_finalize()
  end subroutine rism3d_destroy


  !! PRIVATE



  !> Using the current orientation of the solute, define the minimum box size
  !! and resize all associated grids.
  !!
  !! Periodic : the unit cell is used to determine the number of grid points.
  !!
  !! @param[in] this rism3d object
  subroutine resizeBox(this,ng3)
    use constants_rism, only : PI
    use rism_util, only : isprime, lcm, isFactorable, largestPrimeFactor
    implicit none
#if defined(MPI)
    include 'mpif.h'
    integer :: ierr
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: ng3(3)
    integer :: ngr(3)
    integer :: id
    _REAL_ :: boxlen(3)
    integer :: primes(4) = (/2, 3, 5, 7/)
    logical :: cuv_dimension_changed
    integer :: ngr_size

    ! To properly distribute the run, y and z dimension must have a
    ! number of grid points that is a multiple of mpisize.  Thus,
    ! mpisize must be factorizable by small prime numbers.
    if (.not. isFactorable(this%mpisize, primes)) then
       call rism_report_error("(a,10i4)", &
            "Sorry, 3D-RISM requires that the number " &
            //"of processes be a product of ", primes)
    end if

    ! If we have a fixed box size, we simply retain all of the
    ! previously calculated box size values.

    if (this%periodic) then

       boxlen(:) = this%grid%unitCellLengths(:)

       if ( ng3(1) == -1 .and. ng3(2) == -1 .and. ng3(3) == -1 ) then

          ! Ensure gridpoints fit in unit cell perfectly by adjusting
          ! grid spacing. Current approach treats user-specified grid
          ! spacing as a maximum spacing.
          ngr(:) = ceiling(boxlen(:) / this%grid%spacing(:))

       else

          ! get the grid size from the user input:
          ngr(:) = ng3(:)

       end if

       this%grid%spacing(:) = boxlen(:) / ngr(:)

       ! Determine if the number of grid points in each dimension
       ! product only has prime factors 2, 3, 5 or 7. If not
       ! increment the number of points (in that dimension) until
       ! this is true.

       ! Make sure that each dimension is divisible by 2 and that the
       ! y- and z-dimensions are divisible by this%mpisize if
       ! this%mpisize > 2.  This former is done for the sake of FFT
       ! libraries.
       ! dac note: this code does not ensure that the y and z dimesions
       !   will still be even after dividing things up among the MPI
       !   threads....
       ! dac note: why do both y and z need to be multiples of mpisize?
       !   is in not just z?

       ngr(:) = ngr(:) + mod(ngr(:), 2)
       if (this%mpisize > 2) then
          do id = 2, 3
             if (mod(ngr(id), this%mpisize) /= 0) then
                ngr(id) = ngr(id) + lcm(this%mpisize, 2) &
                     - mod(ngr(id), lcm(this%mpisize, 2))
             end if
          end do
       end if

       do id = 1, 3
          do while (largestPrimeFactor(ngr(id)) .gt. 9)
             if (this%mpisize > 1 .and. id > 1) then
                ngr(id) = ngr(id) + lcm(this%mpisize, 2)
             else
                ngr(id) = ngr(id) + 2
             end if
          end do
       end do
       this%grid%spacing = boxlen / ngr

       cuv_dimension_changed = .false.
       do id = 1, 3
          ngr_size = merge(ngr(id), ngr(id) / this%mpisize, id /= 3)
          if (ubound(this%cuv, id) /= ngr_size &
               .or. ubound(this%oldcuv, id) /= ngr_size) then
             cuv_dimension_changed = .true.
             exit
          end if
       end do
       if (.not. associated(this%cuv) &
            .or. .not. associated(this%oldcuv) &
            .or. cuv_dimension_changed) then
          call reallocateBox(this, ngr, this%grid%spacing)
       end if

       if (this%potential%cutoff >= this%grid%inscribedSphereRadius) then
          call rism_report_error('(a,f8.2,a)', 'solvcut must be < ', &
               this%grid%inscribedSphereRadius, &
               ', the largest inscribed sphere radius of the unit cell.')
       end if

    end if

  end subroutine resizeBox

  !> Using the current box size and resize all associated grids and
  !! variables.
  !! @param[in,out] this rism3d object.
  !! @param[in] ngr Number of grid points along each axis.
  !! @param[in] grdspc Grid spacing along each axis.
  subroutine reallocateBox(this, ngr, grdspc)
    use constants_rism, only : pi
    use rism3d_fft_c
    use safemem
    implicit none
    type(rism3d) :: this
    integer, intent(in) :: ngr(3)
    _REAL_, intent(in) :: grdspc(3)
    integer :: i, id, irank
    _REAL_ :: memuse
    _REAL_ :: unitCellDimensions(6)

    call rism3d_fft_setgrid(this%grid, ngr, grdspc, &
        this%solvent%numAtomTypes, this%fft_aligned)

    !
    ! 2) allocation
    !
    ! reallocate arrays that require preservation of their contents

    ! I THINK WE CAN GET RID OF THE MPI HERE
#if defined(MPI)
    this%cuvWRK => safemem_realloc(this%cuvWRK, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3), &
         this%solvent%numAtomTypes, this%NVec, .true., .true.)
    if (rism3d_solute_charged(this%solute)) then
       this%oldcuvChg => safemem_realloc(this%oldcuvChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvChg
    else
       this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvNoChg
    end if
#else
    this%cuvWRK => safemem_realloc(this%cuvWRK, &
         this%grid%localDimsR(1), this%grid%localDimsR(2), this%grid%localDimsR(3), &
         this%solvent%numAtomTypes, this%NVec, .true., .true.)
    if (rism3d_solute_charged(this%solute)) then
       this%oldcuvChg => safemem_realloc(this%oldcuvChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvChg
    else
       this%oldcuvNoChg => safemem_realloc(this%oldcuvNoChg, this%grid%localDimsR(1), this%grid%localDimsR(2), &
            this%grid%localDimsR(3), this%solvent%numAtomTypes, this%ncuvsteps, .true., .true.)
       this%oldcuv => this%oldcuvNoChg
    end if
#endif /*MPI*/
    ! reallocate arrays that do not require preservation of their contents
    this%guv => safemem_realloc(this%guv, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%huv => safemem_realloc(this%huv, this%grid%totalLocalPointsK, this%solvent%numAtomTypes, &
         o_preserve = .false., o_aligned = .true.)
    this%cuvresWRK => safemem_realloc(this%cuvresWRK, this%grid%totalLocalPointsR, this%solvent%numAtomTypes, &
         this%NVec, .false.)
    this%xvva => safemem_realloc(this%xvva, this%grid%waveNumberArraySize * (this%solvent%numAtomTypes)**2, .false.)

    call rism3d_fft_destroy(this%fft)
    call rism3d_fft_new(this%fft, &
         this%fftw_planner, this%fftw_localtrans, this%fft_aligned, &
         this%grid, &
         this%guv, this%huv)

    ! updated pointers
    call mdiis_resize(this%mdiis_o, this%cuvWRK, this%cuvresWRK, &
         this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%nvec)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))


    ! If we have a previous solution with no charges, then this%cuv
    ! will likely have the wrong previous solution, whether or not we
    ! now have charges set.  The first solution in this%oldcuv will
    ! have the last solution with the current charge setting, so we
    ! copy that over.

    if (this%nsolutionNoChg > 0) then
       if (all(ubound(this%oldcuv(:,:,:,:,1)) .eq. ubound(this%cuv))) then
          call dcopy(product(ubound(this%cuv)), this%oldcuv, 1,this%cuv, 1)
       else
          call rism_report_error('(a,l,a,5(i),a,4(i))', &
               'Should not get here. Attempting to copy oldcuv to cuv failed:'//NEW_LINE('A')//&
               'Charged : ',this%solute%charged,NEW_LINE('A')//&
               'ubound(oldcuv) : ', ubound(this%oldcuv), NEW_LINE('A')//&
               'ubound(cuv) : ', ubound(this%cuv))
      end if
    end if
    
    !
    ! 3) the remaining variables are handled by rism_setup_wavevector, interpolateSolventSusceptibility
    !
    call interpolateSolventSusceptibility(this, this%solvent%xvv, this%xvva)
  end subroutine reallocateBox

  !> Interpolate the solvent-solvent susceptibility, solved on the
  !! 1D-RISM grid, to the 3D-RISM grid.
  !! @param[in,out] this rism3d object.
  !! @param[in] xvv 1D-RISM Xvv
  !! @param[out] xvva Interpolated result.
  subroutine interpolateSolventSusceptibility(this, xvv, xvva)
    use rism_util, only : polynomialInterpolation
    implicit none
    type(rism3d), intent(inout) :: this
    _REAL_, intent(in) :: xvv(:, :, :)
    _REAL_, intent(out) :: xvva(:)
    integer :: iwn, igk, igk1, iv1, iv2
    _REAL_ :: err
    !> Maximum number of points to interpolate.
    integer :: maxPointsToInterp
    parameter (maxPointsToInterp = 5)

    ! Checking R-grid size.
    if (this%grid%waveNumbers(this%grid%waveNumberArraySize) > this%solvent%waveNumbers(this%solvent%numRDFpoints)) then
       call rism_report_error('(a,1pe16.8,a,1pe16.8)', &
            'DISTVV: bulk solvent Kmax=', this%solvent%waveNumbers(this%solvent%numRDFpoints), &
            'insufficient for 3D-grid Kmax=', this%grid%waveNumbers(this%grid%waveNumberArraySize))
    else if (maxPointsToInterp > this%solvent%numRDFpoints) then
       call rism_report_error('(a,i7,a,i7)', &
            'DISTVV: bulk solvent grid size Nr=', this%solvent%numRDFpoints, &
            'insufficient for interpolation maxPointsToInterp=', maxPointsToInterp)
    end if

    ! Interpolate Xvv(k) on 1D-RISM grid to 3D-RISM grid.
    ! Interpolation is performed about the 1D-RISM wave number just
    ! larger than the current 3D-RISM wave number, with a point range
    ! of +/- maxPointsToInterp / 2.
    do iwn = 1, this%grid%waveNumberArraySize
       ! Find the smallest 1D-RISM wave number just larger than the
       ! current 3D-RISM wave number.  The range starts at
       ! this point - maxPointsToInterp/2 and goes up to + maxPointsToInterp/2.
       do igk = 1, this%solvent%numRDFpoints - maxPointsToInterp + 1
          igk1 = igk
          if (this%solvent%waveNumbers(igk1 + maxPointsToInterp/2) > this%grid%waveNumbers(iwn)) then
             exit
          end if
       end do
       ! Interpolate from 1D-RISM to 3D-RISM grids about the midpoint
       ! +/- maxPointsToInterp/2.
       do iv2 = 1, this%solvent%numAtomTypes
          do iv1 = 1, this%solvent%numAtomTypes
             call polynomialInterpolation( &
                  this%solvent%waveNumbers(igk1:igk1 + maxPointsToInterp), &
                  xvv(igk1:igk1 + maxPointsToInterp, iv1, iv2), maxPointsToInterp, &
                  this%grid%waveNumbers(iwn), &
                  xvva(iwn + (iv1 - 1) * this%grid%waveNumberArraySize &
                  + (iv2 - 1) * this%grid%waveNumberArraySize * this%solvent%numAtomTypes), &
                  err)
          end do
       end do
    end do
  end subroutine interpolateSolventSusceptibility

  !!!!
  !! Subroutines to find the iterative 3D-RISM solution.
  !!!!

  !> Main driver for the 3D-RISM solver.
  !! Makes an initial guess of the direct correlation function and
  !! then solve the RISM and closure relations until either the
  !! solution converges or the maximum of steps is reached.
  !! @param[in,out] this rism3d object.
  !! @param[in] ksave Save itermediate results every ksave interations
  !!  (0 means no saves).
  !! @param[in] kshow Print parameter for relaxation steps every kshow
  !!  iteration (0 means no saves).
  !! @param[in] maxSteps Maximum number of rism relaxation steps.
  !! @param[in] tolerance Tolerance in.
  subroutine solve3DRISM(this, ksave, kshow, maxSteps, tolerance)
    use mdiis_c
    use rism3d_restart
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    integer, intent(in) :: ksave, kshow, maxSteps
    _REAL_, intent(in) :: tolerance
    character(72) :: cuvsav = 'rism.csv', guvfile
    integer :: guv_local = 77

    integer :: iatv, igx, igy, igz
    logical :: found, converged = .false.
    integer :: ig, igr, iv, iv2, istep
    _REAL_ :: residual = 0

    ! Absolute first time in solve3DRISM.
    logical, save :: first = .true.

    ! MPI rank counter.
    integer :: irank
    ! iostat
    integer :: stat
    integer :: ientry, nentry
    ! MPI error.
    logical :: ierr

    !if (this%verbose >= 1 .and. size(this%closureList) > 1) &
    !     call rism_report_message("|Using "// &
    !     trim(rism3d_closure_type(this%closure))//" closure")

    ! Make initial guess for DCF.
    this%cuvres = 0
#ifdef MPI
#else
    if (this%mpirank == 0) then
       inquire (file = cuvsav, exist = found)
       if (found .and. first .and. ksave /= 0) then
          write(6,*)'| reading saved Cuv file:  ', cuvsav
          call readRestartFile(cuvsav, this%cuv(:, :, :, :), &
               this%grid%totalLocalPointsR, this%solvent%numAtomTypes)
       else
#endif /*MPI*/

          if (this%nsolution == 0 .or. this%ncuvsteps == 0) then
             ! Default initial guess for DCF is zero everywhere.
             this%cuv(:,:,:,:) = 0.d0
          end if
#ifdef MPI
#else
       end if
    end if
#endif /*MPI*/

    ! Solve 3D-RISM for the current solute-solvent system.
    !
    ! This is done by first using the DCF guess above, the calculated
    ! potentials, and the 1D-RISM solvent-solvent susceptibility to
    ! solve the 3D-RISM equation for the TCF. These values are then
    ! used in the bridge function to obtain a new solution for TCF.
    ! The difference between the bridge function TCF and the 3D-RISM
    ! TCF is the residual, which is used to create a new guess for the
    ! DCF.
    !
    ! The process repeats until the residual is minimized to a desired
    ! amount, at which point the solution is considered converged.
    !
    ! MDIIS is used to increase the rate of convergence by combining
    ! knowledge of previous guesses of the DCF with the current TCF
    ! residual to give an improved guess for the DCF.

    call mdiis_reset(this%mdiis_o)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))

    call timer_start(TIME_R1RISM)
    do istep = 1, maxSteps
       !  -----------------------------------------------------------------
       ! One iteration of 3D-RISM and closure relation, advancing w/ mdiis.
       !  -----------------------------------------------------------------
       call single3DRISMsolution(this, residual, converged, tolerance, &
           this%periodic)

       ! Showing selected and last relaxation steps.
       if (kshow /= 0 .and. this%mpirank == 0 .and. this%verbose >= 2) then
          if (converged .or. mod(istep, kshow) == 0 .or. &
               ksave > 0 .and. mod(istep, max(ksave, 1)) == 0) then
             call rism_report_message('(a,i5,5x,a,1pg10.3,5x,a,i3)', &
                  ' Step=', istep, 'Resid=', residual, 'IS=', &
                  getCurrentNVec(this%mdiis_o))
             call rism_report_flush()
          end if
       end if

       ! Exiting relaxation loop on convergence.
       if (converged) exit
    end do

    if (.not. converged) then
       call rism_report_error('(a,i5)', &
          'RXRISM: reached limit # of relaxation steps: ', maxSteps)
    end if
    call timer_stop(TIME_R1RISM)

    first = .false.
    if (this%mpirank == 0 .and. this%verbose >= 1) then
       call rism_report_message('(a,i5,a)', &
           "|RXRISM converged in ", istep, " steps")
    end if
    return
  end subroutine solve3DRISM


  !> One relaxation step for the UV-RISM equation with the HNC closure,
  !! Guv(r) = exp(-Uuv(r) + Tuv(r) - DelHv0) + DelHv0
  !! Cuv(r) = Guv(r) - 1 - Tvv(r)
  !! Huv(k) = Cuv(k) * (Wvv(k) + Density * Hvv(k))
  !! TuvRes(r) = Huv(r) - Guv(r) - 1
  !! @param[in,out] this A rism3d object.
  !! @param[in,out] residual ???
  !! @param[in,out] converged Returns true if the solution has converged.
  !! @param[in] tolerance Target residual tolerance for convergence.
  subroutine single3DRISMsolution(this, residual, converged, tolerance, &
      periodic)

    use rism3d_fft_c
    use constants_rism, only : PI, FOURPI, omp_num_threads
    implicit none
#include "def_time.h"
#if defined(MPI)
    include 'mpif.h'
#endif /*defined(MPI)*/
    type(rism3d), intent(inout) :: this
    logical, intent(inout) :: converged
    _REAL_, intent(inout) :: residual
    _REAL_, intent(in) :: tolerance
    logical, intent(in) :: periodic
    integer :: iis
    _REAL_ :: earg, tuv0, tvvr
    integer :: istep

    integer ::  ig1, iga, iv, iv1, iv2, igx, igy, igz, igk, ig
#ifdef FFW_THREADS
    integer :: nthreads, totthreads
    integer, external :: OMP_get_max_threads, OMP_get_num_threads
    logical, external :: OMP_get_dynamic, OMP_get_nested
#endif
    integer :: ierr, irank

    ! --------------------------------------------------------------
    ! Cuv(r) is loaded into the guv array.
    ! --------------------------------------------------------------

    call timer_start(TIME_RISMFFTB)
#if defined(MPI)
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   igk = igx + (igy-1) * (this%grid%localDimsR(1) + 2) &
                        + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                   this%guv(igk, iv) = this%cuv(igx, igy, igz, iv)
                end do
             ! Zero out extra space.
             igk = this%grid%localDimsR(1) + 1 + (igy-1) * (this%grid%localDimsR(1) + 2) &
                  + (igz-1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
             this%guv(igk:igk + 1, iv) = 0.d0
          end do
       end do
    end do
#else
!$omp parallel do private(iv,igx,igy,igz,ig1)  &
!$omp&        num_threads(this%solvent%numAtomTypes)
    do iv = 1, this%solvent%numAtomTypes
          do igz = 1, this%grid%localDimsR(3)
             do igy = 1, this%grid%localDimsR(2)
                do igx = 1, this%grid%localDimsR(1)
                   ig1 = igx + (igy-1) * this%grid%localDimsR(1) &
                             + (igz-1) * this%grid%localDimsR(2) &
                                       * this%grid%localDimsR(1)
                   this%guv(ig1, iv) = this%cuv(igx, igy, igz, iv)
                end do
             end do
          end do
       ! Zero out extra space.
       this%guv(this%grid%totalLocalPointsR + 1:this%grid%totalLocalPointsK, iv) = 0.d0
    end do
!$omp end parallel do
#endif /*defined(MPI)*/
    call timer_stop(TIME_RISMFFTB)

    ! --------------------------------------------------------------
    ! Cuv(r) FFT->K.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMFFT)
#if defined(MPI)
    call  rism3d_fft_fwd(this%fft, this%guv)
    this%guv(2:this%grid%totalLocalPointsK:2, :) = &
         -this%guv(2:this%grid%totalLocalPointsK:2, :)
#else
    call  rism3d_fft_fwd(this%fft, this%guv)
#endif /*defined(MPI)*/
    call timer_stop(TIME_RISMFFT)

    ! --------------------------------------------------------------
    ! Huv(k) by RISM.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMHUVK)
!$omp parallel do private(iv1,iv2,ig1,iga) num_threads(omp_num_threads)
    do ig1 = 1, this%grid%totalLocalPointsK
       do iv1 = 1, this%solvent%numAtomTypes
          this%huv(ig1, iv1) = 0d0
          iga = this%grid%waveVectorWaveNumberMap((ig1 + 1) / 2)
          do iv2 = 1, this%solvent%numAtomTypes
             this%huv(ig1, iv1) = this%huv(ig1, iv1) &
                  + this%guv(ig1, iv2) * this%xvva(iga + (iv2 - 1) &
                     * this%grid%waveNumberArraySize &
                  + (iv1 - 1) * this%grid%waveNumberArraySize &
                     * this%solvent%numAtomTypes)
          end do
       end do
    end do
!$omp end parallel do
    call timer_stop(TIME_RISMHUVK)

    ! ---------------------------------------------------------------
    ! Remove the background charge effect
    ! ---------------------------------------------------------------
    if (this%mpirank == 0 .and. this%solute%charged) then
       this%huv(:2, :) = this%huv(:2, :) + this%potential%huvk0(:, :)
    endif

    ! --------------------------------------------------------------
    ! Huv(k) FFT->R.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMFFT)
#if defined(MPI)
    this%huv(2:this%grid%totalLocalPointsK:2, :) = &
         -this%huv(2:this%grid%totalLocalPointsK:2, :)
    call  rism3d_fft_bwd(this%fft, this%huv)
#else
    call  rism3d_fft_bwd(this%fft, this%huv)
#endif /*defined(MPI)*/
    call timer_stop(TIME_RISMFFT)

    call timer_start(TIME_CLOSURE)
    ! --------------------------------------------------------------
    ! Solve the closure for the RDF.
    ! --------------------------------------------------------------
    call rism3d_closure_guv(this%closure, this%guv, this%huv, this%cuv)
    call timer_stop(TIME_CLOSURE)

    ! --------------------------------------------------------------
    ! Calculate TCF residual for use in estimating DCF residual.
    ! --------------------------------------------------------------
    call timer_start(TIME_RISMRESID)
#ifdef MPI
    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             do igx = 1, this%grid%localDimsR(1)
                ig1 = igx + (igy - 1) * this%grid%localDimsR(1) + &
                     (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
                igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                     + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
                this%cuvres(ig1, iv) = this%guv(igk, iv) - 1d0 - this%huv(igk, iv)
             end do
          end do
       end do
    end do
#else
!$omp parallel do num_threads(omp_num_threads)
    do ig=1,this%grid%totalLocalPointsR
       this%cuvres(ig,:) = this%guv(ig,:) - 1d0 - this%huv(ig,:)
    end do
!$omp end parallel do
#endif
    call timer_stop(TIME_RISMRESID)

    ! --------------------------------------------------------------
    ! MDIIS
    ! --------------------------------------------------------------
    call timer_start(TIME_MDIIS)
    call mdiis_advance(this%mdiis_o, residual, converged, tolerance)
    this%cuv => this%cuvWRK(:, :, :, :, mdiis_getWorkVector(this%mdiis_o))
    this%cuvres => this%cuvresWRK(:, :, mdiis_getWorkVector(this%mdiis_o))
    call timer_stop(TIME_MDIIS)

  end subroutine single3DRISMsolution

  ! PROPAGATE PREVIOUS SOLUTIONS

  !> Calculates a new initial guess for CUV based on the final solutions
  !! from previous timesteps.  The maximum number of previous time
  !! steps to use is provided by the user in ncuvsteps.  However, if
  !! there are not enough previous timesteps only nsolution previous
  !! timesteps will be used.
  !!
  !! See section 2.3.1 and eqs. 8-13 of doi:10.1021/ct900460m for
  !! details.
  !! @param[in] this rism3d object.
  subroutine guessDCF(this)
    implicit none
    type(rism3d) :: this
    integer :: iv, n
    n = this%grid%totalLocalPointsR * this%solvent%numAtomTypes
    if (this%ncuvsteps >= 5 .and. this%nsolution >= 5) then
       call dscal(n, 5d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -10d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 10d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, -5d0, this%oldcuv(:, :, :, :, 4), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 1d0, this%oldcuv(:, :, :, :, 5), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 4 .and. this%nsolution >= 4) then
       call dscal(n, 4d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -6d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 4d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 4), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 3 .and. this%nsolution >= 3) then
       call dscal(n, 3d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
       call daxpy(n, 1d0, this%oldcuv(:, :, :, :, 3), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps >= 2 .and. this%nsolution >= 2) then
       call dscal(n, 2d0, this%cuv(:, :, :, :), 1)
       call daxpy(n, -1d0, this%oldcuv(:, :, :, :, 2), 1, this%cuv(:, :, :, :), 1)
    else if (this%ncuvsteps == 0) then
       this%cuv(:, :, :, :) = 0
    end if
  end subroutine guessDCF


  !> Updates the values in the this%oldcuv queue.  The oldest value (the
  !! ncuvstep index) is pushed out, the remainder of the data is
  !! shifted and the newest solution is placed in the first index.
  !! @param[in] this rism3d object.
  subroutine updateDCFguessHistory(this)
    implicit none
    type(rism3d) :: this
    integer :: istep, iv
#ifdef RISM_DEBUG
    write(0, *) "DCF_GUESS_HISTORY_UPDATE"; call flush(0)
#endif /*RISM_DEBUG*/
    if (this%ncuvsteps == 0) return
    do istep = min(this%ncuvsteps, this%nsolution), 2, -1
       call dcopy(this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%oldcuv(:, :, :, :, istep-1), 1, &
            this%oldcuv(:, :, :, :, istep), 1)
    end do
    call dcopy(this%grid%totalLocalPointsR * this%solvent%numAtomTypes, this%cuv(:, :, :, :), 1, &
         this%oldcuv(:, :, :, :, 1), 1)
  end subroutine updateDCFguessHistory

  subroutine check_xvv_info( this )
    use constants_rism, only : PI, FOURPI
    implicit none
    type(rism3d_potential), intent(inout) :: this !< potential object.

    _REAL_ soluteQ
    logical :: first = .true.

    ! Allocate huvk0 if not done already.
    if (.not. associated(this%huvk0)) then
       this%huvk0 => safemem_realloc(this%huvk0, 2, this%solvent%numAtomTypes)
    end if

    ! Check if the background charge correction is provided. Old
    ! Xvv files only have delhv0, which combines the background
    ! correction and the long range asymptotics. If we only have
    ! delhv0, remove the long-range asymptotics contribution.
    soluteQ = sum(this%solute%charge)/this%grid%boxVolume
    if(all(this%solvent%background_correction .ne. HUGE(1d0))) then
       this%huvk0(1, :) = this%solvent%background_correction(:) * soluteQ
#if 0
       if (first .and. this%grid%offsetK(3) == 0) then
          write(6,'(a,6f10.5)') '|  huvk0 = ', this%huvk0(1, :)
          first = .false.
       end if
#endif
    else
       ! for pure water, kappa is zero, and there should be no
       !    background correction:
       if( this%solvent%xappa == 0.d0 ) then
          this%huvk0(1, :) = 0.d0
       else
          this%huvk0(1, :) = ( this%solvent%delhv0(:) &
             - this%solvent%charge_sp(:) / this%solvent%dielconst &
             * FOURPI / this%solvent%xappa**2 ) * soluteQ 
       end if
    end if
    this%huvk0(2, :) = 0d0
  end subroutine check_xvv_info

end module rism3d_c
