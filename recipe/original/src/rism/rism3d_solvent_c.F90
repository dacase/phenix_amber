!<compile=optimized>

#include "../include/dprec.fh"

!> Object for solvent information coming from 1D-RISM to be used in
!! 3D-RISM.
module rism3d_solvent_c
  use rism_report_c
  use safemem
  implicit none
  type rism3d_solvent     
     !> Solvent temperature. [K]
     _REAL_ :: temperature = 0
     !> Solvent dielectric constant.
     _REAL_ :: dielconst = 0
     !> Inverse Debye length. [1/A]
     _REAL_ :: xappa = 0
     !> Compressibility. [A^3]
     _REAL_ :: xikt = 0
     !> Charge smear parameter required for long range
     !! electrostatics. [A]
     _REAL_ :: smear = 0
     !> Number of solvent atom types.
     integer :: numAtomTypes = 0
     !> Number of molecular species.
     integer :: numMolecules = 0
     !> The number of points in the solvent-solvent radial
     !! distribution function.
     integer :: numRDFpoints = 0

     !> The multiplicity of each atom in the solvent model.
     integer, pointer :: atomMultiplicity(:) => NULL()
     !> Number of solvent atoms (sites) for each molecular species.
     integer, pointer :: numAtoms(:) => NULL()

     !> Names of solvent atoms. Used for the solute-solvent RDF, DCF
     !! and TCF output files.
     character(len=4), pointer :: atomName(:) => NULL()

     !> Grid spacing from 1D-RISM calculation. [A]
     _REAL_ :: gridSpacingR = 0
     !> Fourier grid spacing from 1D-RISM calculation,
     !! dt = pi / (nr * dr). [1/A]
     !! Currently this is used only to set up the solvent
     !! calculations.
     _REAL_ :: gridSpacingK = 0
     !> Precalculated table of the radial grid point positions for
     !! 1D-RISM in Fourier space. (Basically the distance in Fourier
     !! space each grid point is from the origin.)
     !! This is then used for interpolating Xvv for the needs of
     !! 3D-RISM.
     _REAL_, pointer :: waveNumbers(:) => NULL()
     
     !> Solvent-solvent susceptibility.
     _REAL_, pointer :: xvv(:,:,:) => NULL()
     !> Solvent atom charge by atom type. [sqrt(kT A)]
     _REAL_, pointer :: charge(:) => NULL()
     !> Charge of the parent molecular species for each type. [sqrt(kT A)]
     _REAL_, pointer :: charge_sp(:) => NULL()
     !> Solvent density by atom type. [#/A^3]
     _REAL_, pointer :: density(:) => NULL()
     !> Solvent density by molecular species. [#/A^3]
     _REAL_, pointer :: density_sp(:) => NULL()
     !> Solvent Lennard-Jones sigma coefficient by atom type. [A]
     _REAL_, pointer :: ljSigma(:) => NULL()
     !> Solvent Lennard-Jones epsilon coefficient by atom type. [kT]
     _REAL_, pointer :: ljEpsilon(:) => NULL()
     !> coordinates.  Labeled by dimension, multiplicity number, atom number, species
     _REAL_, pointer :: coord(:,:,:,:) => NULL()
     !> Contribution to the TCF at k=0 for both periodic and
     !non-periodic calculations: ! See eq. 20 in Kolavenko/Hirata
     !2000.
     !! -Lim_k->0 ( Sum_v1 Qv1 * Xv1v2(k) 4pi / k^2 ) from 1D-RISM.
     !! This is the same as delhv0 except that it does not include the
     !! TCF long-range asymptotics
     _REAL_, pointer :: background_correction(:) => NULL()
     !> The contribution to the TCF by the background charge which
     !! gives the supercell electroneutrality.  This portion of the term
     !! is only the contribution by the solvent, and does not include the
     !! solute specific terms (box volume and total solute change).
     !! See eq. 20 in Kolavenko/Hirata 2000.
     !! -Lim_k->0 ( Sum_v1 Qv1 * Xv1v2(k) 4pi / k^2 - hlkv0 ) from 1D-RISM.
     _REAL_, pointer :: delhv0(:) => NULL()

     !> True if the any of the solvent species has non-zero charge.
     logical :: ionic

     !> Version number of 1D-RISM solvent-solvent susceptibility (XVV)
     !! file.
     real :: xvv_version
  end type rism3d_solvent

  interface rism3d_solvent_new
     module procedure rism3d_solvent_new_all, rism3d_solvent_new_readxvv
  end interface rism3d_solvent_new

  public :: rism3d_solvent_new, rism3d_solvent_clone, rism3d_solvent_destroy
#ifdef MPI
  public :: rism3d_solvent_mpi_clone
#endif /*MPI*/
  
  private :: allocate_solv, calculateWavenumbers, readxvv1, readxvv2

contains
  
  !> Constructor.  Creates new solvent object be specifying all
  !! data. if this is an MPI run, only the rank 0 process needs valid
  !! arguments (other than the communicator).
  !! @param[in] this The new solvent object.
  !! @param[in] dr Grid spacing from RISM1D calculation [A].
  !! @param[in] dt Fourier grid spacing from RISM1D calculation
  !!     dt = pi/nr/dr [1/A].  Currently this is used only to set up the
  !!     solvent calculations.
  !! @param[in] temperatur Solvent temperature [K].
  !! @param[in] dielconst Solvent dielectric constant.
  !! @param[in] xappa Inverse Debye length [1/A].
  !! @param[in] xikt Compressibility [A^3].
  !! @param[in] smear ???
  !! @param[in] natom Number of solvent atom types.
  !! @param[in] nspecies Number of molecular species.
  !! @param[in] numAtomSpecies Number of solvent sites for each molecular species.
  !! @param[in] nr The number of points in the solvent-solvent RDF.
  !! @param[in] mult The multiplicity of each atom in the solvent model.
  !! @param[in] atomName Name of the solvent atoms.  Used for the GUV,CUV and HUV output files.
  !! @param[in] xvv Original solvent chi.
  !! @param[in] fourier_tbl Precalculated table of the radial spacing
  !!     for 1D-RISM in Fourier space.  dt = pi / (nr * dr) This is then
  !!     used for interpolating Xvv for the needs of 3D-RISM.
  !! @param[in] charge Solvent atom charge by type [sqrt(kT A)].
  !! @param[in] charge_sp Total solvent species charge for each atom by type [sqrt(kT A)].
  !! @param[in] density Solvent density by atom type [#/A^3].
  !! @param[in] density_sp Solvent density by molecular species [#/A^3].
  !! @param[in] ljSigma Solvent LJ-sigma by atom type [A].
  !! @param[in] eps Solvent LJ-epsilon by atom type [kT].
  !! @param[in] background_correction Background charge correction of the solute-solvent TCF
  !!    at k = 0 (not yet multiplied by solute-specific coefficients - see
  !!    eq. 20 in Kovalenko 2000).
  !! @param[in] delhv0 Long range asymptotics and background charge correction of the solute-solvent TCF
  !!    at k = 0 (not yet multiplied by solute-specific coefficients - see
  !!    eq. 20 in Kovalenko 2000).
  !! @param[in] o_mpicomm (optional) MPI communicator.
  subroutine rism3d_solvent_new_all(this, dr, dt, temperature, dielconst, xappa, &
       xikt, smear, natom, nspecies, numAtomSpecies, nr, mult, atomName, &
       fourier_tbl, xvv, charge, charge_sp, density, density_sp, ljSigma, eps, &
       background_correction, delhv0, coord, xvv_version, o_mpicomm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_solvent), intent(inout) :: this
    _REAL_, intent(in) :: dr, dt, temperature, dielconst, xappa, xikt, smear
    integer, intent(in) ::  natom, nspecies, nr
    integer, intent(in) :: mult(natom), numAtomSpecies(nspecies)
    character(len = 4), intent(in) :: atomName(natom)
    _REAL_, intent(in) :: fourier_tbl(nr), xvv(nr, natom, natom), &
         charge(natom), charge_sp(natom), density(natom), density_sp(nspecies), &
         ljSigma(natom), eps(natom), background_correction(natom), delhv0(natom), &
         coord(3,maxval(mult),maxval(numAtomSpecies),nspecies)
    integer, optional, intent(in) :: o_mpicomm
    real, intent(in) :: xvv_version
    integer :: mpicomm, mpirank, mpisize, err


    mpicomm =0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if (present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if (mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLVENT: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if (err /=0) call rism_report_error&
            ("(a,i8)","RISM3D INIT: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if (err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLV: could not get MPI size for communicator ",mpisize)
    end if
#endif /*MPI*/
    if (mpirank == 0) then
       this%xvv_version = xvv_version
       this%gridSpacingR = dr
       this%gridSpacingK = dt
       this%temperature = temperature
       this%dielconst = dielconst
       this%xappa = xappa
       this%xikt = xikt
       this%smear = smear
       this%numAtomTypes = natom
       this%numMolecules = nspecies
       this%numRDFpoints = nr
       call allocate_solv(this)

       ! After allocate so it is not overwritten with HUGE.
       this%atomMultiplicity = mult
       this%numAtoms = numAtomSpecies
       this%atomName = atomName
       this%waveNumbers = fourier_tbl
       this%xvv = xvv
       this%charge = charge
       this%charge_sp = charge_sp
       this%density = density
       this%density_sp = density_sp
       this%ljSigma = ljSigma
       this%ljEpsilon = eps
       this%background_correction = background_correction
       this%delhv0 = delhv0
       
       if (this%xvv_version >= 1.) then
          ! coordinates need to be allocated separately
          this%coord => safemem_realloc(this%coord, 3, maxval(this%atomMultiplicity), maxval(this%numAtoms), &
               this%numMolecules, .true., .false.)
          this%coord = coord
       end if
       call calculateWavenumbers(this)
    end if
#ifdef MPI
    if (mpisize /= 1)&
         call rism3d_solvent_mpi_clone(this, mpirank,mpicomm)
#endif /*MPI*/
    if (sum(abs(this%charge_sp)) > 0d0) then
       this%ionic = .true.
    else
       this%ionic = .false.
    end if
  end subroutine rism3d_solvent_new_all

  
  !> Constructor.  Creates new solvent object directly from a legacy
  !! Amber7 format XVV file. If this is an MPI run, only the rank 0
  !! process needs a valid Xvv file.
  !! @param[in,out] this The new solvent object.
  !! @param[in] xvvfile Name of the XVV file to read.
  !! @param[in] o_mpicomm (optional) MPI communicator.
  subroutine rism3d_solvent_new_readxvv(this, xvvfile, o_mpicomm)
    use rism_util, only : freeUnit
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_solvent),intent(inout) :: this
    character(*), intent(in) :: xvvfile
    integer, optional, intent(in) :: o_mpicomm
    integer :: mpicomm, mpirank, mpisize, err, unit
    integer :: iatom
    mpicomm = 0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if (present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if (mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLVENT: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if (err /= 0) call rism_report_error&
            ("(a,i8)","RISM3D INIT: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if (err /= 0) call rism_report_error&
            ("(a,i8)","RISM3D SOLV: could not get MPI size for communicator ",mpisize)
    end if
#endif /*MPI*/
    if (mpirank == 0) then

       unit = freeUnit()
       !unit conversion done in readxvv()
       open (unit,file=xvvfile,status='old')
       call readxvv1(this,unit,this%numAtomTypes,this%numMolecules,this%numRDFpoints)
       call allocate_solv(this)
       call readxvv2(this,unit)
       close(unit)
       call calculateWavenumbers(this)
       !remove numerical error from reading file
       do iatom = 1, this%numAtomTypes
          if (abs(this%charge_sp(iatom)) < 1d-6) &
               this%charge_sp(iatom) = 0
       end do

#ifdef MPI
       if (mpisize /= 1) &
            call rism3d_solvent_mpi_clone(this, mpirank,mpicomm)
#endif /*MPI*/
    end if

    if (sum(abs(this%charge_sp)) > 0d0) then
       this%ionic = .true.
    else
       this%ionic = .false.
    end if
  end subroutine rism3d_solvent_new_readxvv


!! Clone constructor.  Creates new solvent object that is identical to the
!! original.
!! IN:
!!    this :: object to be copied
!!    clone :: clone of the object.  No memory space is shared.
  subroutine rism3d_solvent_clone(this,clone)
    implicit none
    type(rism3d_solvent),intent(in) :: this
    type(rism3d_solvent),intent(inout) :: clone
    call rism3d_solvent_new(clone, this%gridSpacingR, this%gridSpacingK, this%temperature, this%dielconst, &
         this%xappa, this%xikt, this%smear, &
         this%numAtomTypes, this%numMolecules, this%numAtoms, this%numRDFpoints, this%atomMultiplicity, &
         this%atomName, this%waveNumbers, this%xvv, this%charge, this%charge_sp, &
         this%density, this%density_sp, &
         this%ljSigma, this%ljEpsilon, this%background_correction, this%delhv0, &
         this%coord, this%xvv_version)
  end subroutine rism3d_solvent_clone

#ifdef MPI

!! Allocates memory on non-master nodes and then distributes information out
!! from the master.  It is assumed that the object on the master already exists.

  subroutine rism3d_solvent_mpi_clone(this,rank,comm)
    implicit none
    type(rism3d_solvent),intent(inout) :: this
    integer, intent(in) :: rank,comm
    integer :: err
    include 'mpif.h'
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%xvv_version, 1, mpi_real, 0, comm, err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast XVV_VERSION")
    call mpi_bcast(this%numRDFpoints,1,mpi_integer,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast NR")
    call mpi_bcast(this%numAtomTypes,1,mpi_integer,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast NATOM")
    call mpi_bcast(this%numMolecules,1,mpi_integer,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast NSPECIES")

    ! non-master processes should now allocate memory
    ! see below for COORD, which requires MULTIPLICITY to be broadcast first
    if (rank /= 0) then
       call allocate_solv(this)
    end if

    !now distribute the arrays to the non-master processes
    call mpi_bcast(this%temperature,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast TEMPERATURE")
    call mpi_bcast(this%dielconst,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast DIELCONST")
    call mpi_bcast(this%xappa,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast XAPPA")
    call mpi_bcast(this%xikt,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast XIKT")
    call mpi_bcast(this%smear,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast SMEAR")
    call mpi_bcast(this%gridSpacingR,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast GRIDSPACINGR")
    call mpi_bcast(this%gridSpacingK,1,mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast GRIDSPACINGK")
    call mpi_bcast(this%xvv,size(this%xvv,1)*size(this%xvv,3)&
         *size(this%xvv,3),mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast XVV")
    call mpi_bcast(this%waveNumbers,size(this%waveNumbers),&
         mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast FOURIER_TBL")
    call mpi_bcast(this%charge,size(this%charge),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast CHARGE")
    call mpi_bcast(this%charge_sp,size(this%charge_sp),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast CHARGE_SP")
    call mpi_bcast(this%background_correction,size(this%background_correction),mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast BACKGROUND_CORRECTION")
    call mpi_bcast(this%delhv0,size(this%delhv0),mpi_double_precision,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast DELHV0")
    call mpi_bcast(this%density,size(this%density),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast DENSITY")
    call mpi_bcast(this%density_sp,size(this%density_sp),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DENSITY_SP")
    call mpi_bcast(this%ljSigma,size(this%ljSigma),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast SIG_S")
    call mpi_bcast(this%ljEpsilon,size(this%ljEpsilon),mpi_double_precision,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast EPS")

    call mpi_bcast(this%atomMultiplicity,size(this%atomMultiplicity),mpi_integer,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast MULT")
    call mpi_bcast(this%numAtoms,size(this%numAtoms),mpi_integer,0,&
         comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast NATOMSPECIES")
    call mpi_bcast(this%atomName,size(this%atomName)*4,mpi_character,0,comm,err)
    if (err /=0) call rism_report_error&
         ("RISM3D_SOLVENT: could not broadcast ATOMNAME")

    ! allocate COORD memory here
    if(this%xvv_version >= 1.00) then
       if (rank /= 0) then
          this%coord => safemem_realloc(this%coord, 3, maxval(this%atomMultiplicity), maxval(this%numAtoms), &
               this%numMolecules, .true., .false.)
       end if
       call mpi_bcast(this%coord,size(this%coord),mpi_double_precision,0,comm,err)
       if (err /=0) call rism_report_error&
            ("RISM3D_SOLVENT: could not broadcast COORD")
    end if
    
    !non-master nodes finish initialization
    if (rank /=0) then
       call calculateWavenumbers(this)
    end if

    !this is cheaper to recalculate than to broadcast
    if (sum(abs(this%charge_sp)) > 0d0) then
       this%ionic = .true.
    else
       this%ionic=.false.
    end if
  end subroutine rism3d_solvent_mpi_clone
#endif /*MPI*/

!! destroyer

  subroutine rism3d_solvent_destroy(this)
    use safemem
    implicit none
    type(rism3d_solvent) :: this
    integer :: err
    if (safemem_dealloc(this%atomMultiplicity)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate MULT")
    if (safemem_dealloc(this%numAtoms)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate NATOMSPECIES")
    if (safemem_dealloc(this%atomName)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate ATOMNAME")
    if (safemem_dealloc(this%waveNumbers)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate WAVENUMBERS")
    if (safemem_dealloc(this%xvv)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate XVV")
    if (safemem_dealloc(this%charge)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate CHARGE")
    if (safemem_dealloc(this%charge_sp)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate CHARGE_SP")
    if (safemem_dealloc(this%density)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate DENSITY")
    if (safemem_dealloc(this%density_sp)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate DENSITY_SP")
    if (safemem_dealloc(this%ljSigma)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate LJSIGMA")
    if (safemem_dealloc(this%ljEpsilon)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate LJEPSILON")
    if (safemem_dealloc(this%background_correction)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate BACKGROUND_CORRECTION")
    if (safemem_dealloc(this%delhv0)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate DELHV0")
    if (safemem_dealloc(this%coord)/=0)&
         call rism_report_error("RISM3D_SOLVENT: failed to deallocate COORD")
  end subroutine rism3d_solvent_destroy


!!  PRIVATE



  !> Allocate memory at the begining of the simulation.
  !! @param[in,out] this solvent object.
  subroutine allocate_solv(this)
    use safemem
    implicit none

    type(rism3d_solvent), intent(inout) :: this
    ! New stack limit.
    integer :: new_limit    

    ! Allocate real memory.
    this%waveNumbers => safemem_realloc(this%waveNumbers,this%numRDFpoints,.false.)
    this%xvv => safemem_realloc(this%xvv,this%numRDFpoints,this%numAtomTypes,this%numAtomTypes,.false.)
    this%charge => safemem_realloc(this%charge,this%numAtomTypes,.false.)
    this%charge_sp => safemem_realloc(this%charge_sp,this%numAtomTypes,.false.)
    this%ljEpsilon => safemem_realloc(this%ljEpsilon,this%numAtomTypes,.false.)
    this%ljSigma => safemem_realloc(this%ljSigma,this%numAtomTypes,.false.)
    this%density => safemem_realloc(this%density,this%numAtomTypes,.false.)
    this%density_sp => safemem_realloc(this%density_sp,this%numMolecules,.false.)
    this%background_correction => safemem_realloc(this%background_correction,this%numAtomTypes,.false.)
    this%delhv0 => safemem_realloc(this%delhv0,this%numAtomTypes,.false.)

    ! Allocate integer memory.
    this%atomMultiplicity => safemem_realloc(this%atomMultiplicity,this%numAtomTypes,.false.)
    this%numAtoms => safemem_realloc(this%numAtoms, this%numMolecules, .false.)

    ! Allocate character (holrith) memory.
    this%atomName => safemem_realloc(this%atomName,4, this%numAtomTypes,.false.)

  end subroutine allocate_solv


  !> Pre-calculate reciprocal grid spacing and wave numbers.
  !! @param[in,out] this solvent object.
  subroutine calculateWavenumbers(this)
    use constants_rism, only : PI
    implicit none
    type(rism3d_solvent), intent(inout) :: this
    integer :: ir
#ifdef RISM_DEBUG
    write(6,*) "CALCULATEWAVENUMBERS"; call flush(6)
#endif /*RISM_DEBUG*/
    ! Compute the 1D Fourier grid spacing used for the 1D-RISM calculations.
    this%gridSpacingK = PI / (this%numRDFpoints * this%gridSpacingR)
    do ir = 1, this%numRDFpoints
       this%waveNumbers(ir) = this%gridSpacingK * (ir - 1)
    end do
  end subroutine calculateWavenumbers


!! Reads in the header data from a Amber style Xvv file.  Determines the version
!! and reads in POINTER information that determines the sizes of array in the
!! rest of the file
!! IN:
!!    nf       : Fortran unit number for an open file
!!    natom    : total number of solvent sites
!!    nspecies : number of molecular species
!!    nr       : number of grid points in Xvv

  subroutine readxvv1(this,nf,natom,nspecies,nr)
    use rism_parm
    implicit none
    type(rism3d_solvent), intent(inout) :: this
    integer, intent(in) :: nf
    integer, intent(out) :: natom, nspecies, nr
    integer :: i, iok
    character(len=80) fmt, filename
    character(len=80) fmtin,ifmt,afmt,rfmt,type

    inquire(unit=nf,name=filename)
    fmt = ''
    ifmt = '(10I8)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    call rism_parm_nxtsec_init()
    !     ----- READ THE POINTERS AND THE VERSION -----

    fmtin = ifmt
    type = 'POINTERS'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok,this%xvv_version)
    !check for a version flag
    if (iok==-1) call rism_report_error('Xvv file '//trim(filename)//&
         ': missing %VERSION line')
    !check the version number.  Acceptable values are 0000.001, 0001.000 and 0001.001
    if (this%xvv_version == 0.000d0 .or. &
         (this%xvv_version >= 0.0015d0 .and. this%xvv_version <= 0.9995) .or. &
         this%xvv_version > 1.0015d0) then
       call rism_report_error('Xvv file '//trim(filename)//&
            ': bad version number. Must be 0001.001, 0001.000 or 0000.001')
    else if (this%xvv_version < 1.0005d0) then
       call rism_report_warn("(a,f7.3)","Xvv version: ",dble(this%xvv_version))
       call rism_report_warn("Unable to calculate UC or PMV temperature derivatives")
       if (this%xvv_version < 0.9995d0) then
          call rism_report_warn("Unable to calculate energy/entropy decomposition")
          call rism_report_warn("UC assumes pure water")
       end if
    end if

    read(nf,fmt) NR,NATOM,NSPECIES

  end subroutine readxvv1


!! Reads the data from a Amber style Xvv file into a solvent object.  Does
!! necessary unit conversion. Number of sites and grid size should be set in the
!! solvent object based on previous call to readxvv1().
!! IN:
!!    solv : solvent object
!!    nf   : Fortran unit number of an open file

  subroutine readxvv2(this, nf)
    use constants_rism, only : COULOMB_CONST_E, KB, BOLTZMANN, AVOGADRO
    use rism_parm
    implicit none
    type(rism3d_solvent), intent(inout) :: this
    integer, intent(in) :: nf

    integer :: i, iok
    character(len=80) fmt, filename
    character(len=80) fmtin, ifmt, afmt, rfmt, type

    integer itab, itab1, ig, iv, imv, iv1, iv2, iiga, isp


    inquire(unit=nf, name=filename)
    fmt = ''

    ifmt = '(10I8)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----

    fmtin = rfmt
    type = 'THERMO'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) this%temperature,this%dielconst,this%xappa,this%xikt,this%gridSpacingR,this%smear

    fmtin = ifmt
    type = 'MTV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%atomMultiplicity(i),i = 1,this%numAtomTypes)

    if (this%xvv_version >= 1.) then
       fmtin = ifmt
       type = 'NVSP'
       call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
       read(nf, fmt) (this%numAtoms(i), i = 1, this%numMolecules)
    else
       this%numAtoms=huge(1)
    end if

    fmtin = afmt
    type = 'ATOM_NAME'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%atomName(i),i = 1,this%numAtomTypes)

    fmtin = rfmt
    type = 'RHOV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%density(i),i = 1,this%numAtomTypes)

    fmtin = rfmt
    type = 'QSPV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%charge_sp(i),i = 1,this%numAtomTypes)
    ! UNIT CONVERSION
    ! from [e] to [sqrt(kT A)].
    if (this%xvv_version < 1.) then
       this%charge_sp = this%charge_sp &
            * sqrt(COULOMB_CONST_E / (KB *this%temperature))
    end if

    fmtin = rfmt
    type = 'QV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0, fmtin,  type,  fmt,  iok)
    read(nf, fmt) (this%charge(i), i = 1, this%numAtomTypes)
    ! UNIT CONVERSION
    ! from [e] to [sqrt(kT A)].
    if (this%xvv_version < 1.) then
       this%charge = this%charge &
            * sqrt(COULOMB_CONST_E / (KB * this%temperature))
    end if

    fmtin = rfmt
    type = 'EPSV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%ljEpsilon(i),i = 1,this%numAtomTypes)
    ! UNIT CONVERSION
    ! from [J/MOL] to [kT].
    !  this%ljEpsilon=this%ljEpsilon/JPKC
    if (this%xvv_version < 1.) then
       this%ljEpsilon = this%ljEpsilon / (BOLTZMANN * AVOGADRO * this%temperature)
    end if

    fmtin = rfmt
    if (this%xvv_version < 1.) then
       type = 'SIGV'
    else
       type = 'RMIN2V'
    end if
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%ljSigma(i),i = 1,this%numAtomTypes)

    ! this is optional
    fmtin = rfmt
    type = 'BACKGROUND_CORRECTION'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 1, fmtin,  type,  fmt,  iok)
    if (iok == 0) then
       read(nf,fmt) (this%background_correction(i),i = 1,this%numAtomTypes)
    else
       this%background_correction = HUGE(1d0)
    end if

    fmtin = rfmt
    type = 'DELHV0'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%delhv0(i),i = 1,this%numAtomTypes)
    !UNIT CONVERSION
    !from [e] to [sqrt(kT A)]
    if (this%xvv_version<1.)&
         this%delhv0 = this%delhv0/sqrt(COULOMB_CONST_E/KB/this%temperature)


    if (this%xvv_version >= 1.) then
       this%coord => safemem_realloc(this%coord, 3, maxval(this%atomMultiplicity), maxval(this%numAtoms), &
            this%numMolecules, .true., .false.)
       fmtin = rfmt
       type = 'COORD'
       call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
       iv2=0
       do isp = 1, this%numMolecules
          do iv1 = 1, this%numAtoms(isp)
             iv2 = iv2+1
             do imv = 1, this%atomMultiplicity(iv2)
                read(nf,fmt) (this%coord(i,imv,iv1,isp),i = 1,3)
             end do
          end do
       end do
    end if
    
    fmtin = rfmt
    type = 'XVV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (((this%xvv(itab,iv1,iv2),itab = 1,this%numRDFpoints),iv1=1,this%numAtomTypes),iv2=1,this%numAtomTypes)

    ! Fields for version 1.000.
    if (this%xvv_version >= 1.) then
       ! Get DENSITYSP from other information.
       do isp = 1, this%numMolecules
          this%density_sp(isp) = this%density(sum(this%numAtoms(1:isp))) / &
               this%atomMultiplicity(sum(this%numAtoms(1:isp)))
       end do
    else
       this%density_sp = 0d0
       this%density_sp(1) = this%density(1)
    end if
  end subroutine readxvv2
end module rism3d_solvent_c
