! Do not optimize right now -- segfaults with ifort

! This file contains the code necessary to initialize the sander state without
! reading command-line flags or input files, as well as computing energies and
! forces via an API.
#include "../include/assert.fh"

#ifdef LES
module sanderles_api
#else
module sander_api
#endif

   use constants, only: MAX_QUANTUM_ATOMS
   use file_io_dat, only: MAX_FN_LEN
   use prmtop_type, only: prmtop_struct, destroy_prmtop_struct, read_prmtop_file
   use qmmm_module, only: qmmm_input_options
   use state, only: potential_energy_rec

! qmmm_input options is a data structure defined in qmmm_module, shown below:
!
! type qmmm_input_options
!    ! Allow a way to input options programmatically through an API rather than
!    ! requiring an input file
!    sequence
!    integer :: iqmatoms(MAX_QUANTUM_ATOMS), qmgb, lnk_atomic_no, &
!               ndiis_matrices, ndiis_attempts, lnk_method, qmcharge, &
!               corecharge, buffercharge, spin, qmqmdx, verbosity, &
!               printcharges, printdipole, print_eigenvalues, peptide_corr, &
!               itrmax, printbondorders, qmshake, qmmmrij_incore, &
!               qmqm_erep_incore, pseudo_diag, qm_ewald, qm_pme, kmaxqx, &
!               kmaxqy, kmaxqz, ksqmaxq, kappa, qmmm_int, adjust_q, tight_p_conv, &
!               diag_routine, density_predict, fock_predict, vsolv, &
!               dftb_maxiter, dftb_disper, dftb_chg, abfqmmm, hot_spot, &
!               qmmm_switch, core_iqmatoms(MAX_QUANTUM_ATOMS), &
!               buffer_iqmatoms(MAX_QUANTUM_ATOMS)
!    _REAL_ :: qmcut, lnk_dis, scfconv, errconv, dftb_telec, dftb_telec_step, &
!              fockp_d1, fockp_d2, fockp_d3, fockp_d4, damp, vshift, &
!              pseudo_diag_criteria, min_heavy_mass, r_switch_hi, r_switch_lo
!    character(len=8192) :: qmmask, coremask, buffermask, centermask
!    character(len=256) :: dftb_3rd_order, dftb_slko_path
!    character(len=12) :: qm_theory
! end type qmmm_input_options

! potential_energy_rec is a data structure defined in state.F90, shown below:

! type potential_energy_rec
!   sequence
!   _REAL_  :: tot
!   _REAL_  :: vdw
!   _REAL_  :: elec
!   _REAL_  :: gb
!   _REAL_  :: bond
!   _REAL_  :: angle
!   _REAL_  :: dihedral
!   _REAL_  :: vdw_14
!   _REAL_  :: elec_14
!   _REAL_  :: constraint
!   _REAL_  :: polar
!   _REAL_  :: hbond
!   _REAL_  :: surf
!   _REAL_  :: scf
!   _REAL_  :: disp
!   _REAL_  :: dvdl
!   _REAL_  :: angle_ub
!   _REAL_  :: imp
!   _REAL_  :: cmap
!   _REAL_  :: emap
!   _REAL_  :: les
!   _REAL_  :: noe
!   _REAL_  :: pb
!   _REAL_  :: rism
!   _REAL_  :: ct
!   _REAL_  :: amd_boost
! end type potential_energy_rec

   ! Datatypes
   type :: sander_input
      sequence
      ! Floats (input parameters)
      double precision :: extdiel
      double precision :: intdiel
      double precision :: rgbmax
      double precision :: saltcon
      double precision :: cut
      double precision :: dielc
      double precision :: rdt
      double precision :: fswitch
      double precision :: restraint_wt
      double precision :: grdspc1
      double precision :: mdiis_del

      ! Integers (toggle options)
      integer :: igb
      integer :: alpb
      integer :: gbsa
      integer :: lj1264
      integer :: ipb
      integer :: inp
      integer :: vdwmeth
      integer :: ew_type
      integer :: ntb
      integer :: ifqnt
      integer :: irism
      integer :: rism_verbose
      integer :: jfastw
      integer :: ntf
      integer :: ntc
      integer :: ntr
      integer :: ibelly
      integer :: mask_from_ref ! Antoine Marion

      ! Character string arrays
      character(len=256) :: restraintmask
      character(len=256) :: bellymask
      character(len=256) :: refc
   end type sander_input

   private

   logical, save :: is_setup_ = .false.

   integer, dimension(:), allocatable :: ipairs

   public :: gas_sander_input, pme_sander_input, sander_setup, energy_forces, &
             set_positions, set_box, sander_input, qmmm_input_options, &
             potential_energy_rec, sander_cleanup, MAX_FN_LEN, get_box, &
             qm_sander_input, sander_natom, prmtop_struct, read_inpcrd_file, &
             get_inpcrd_natom, destroy_prmtop_struct, read_prmtop_file, &
             is_setup, get_positions

contains

! Returns .true. if sander is set up and .false. otherwise
logical function is_setup()

   implicit none

   is_setup = is_setup_

   return

end function is_setup

! Sets up input options for a gas-phase (i.e., no periodic boundaries) system
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
! gb : integer (optional)
!     If GB model is desired, igb will be set to this.  If not present or 0, a
!     gas phase calculation is set up. If set to 10, a PB calculation is set up
!
! If gb is given, the given GB model is set up. Otherwise, the system is set up
! for a gas-phase simulation
subroutine gas_sander_input(inp, gb)

   use constants, only: NO_INPUT_VALUE

   implicit none

   type(sander_input), intent(out) :: inp
   integer, optional, intent(in)   :: gb

   inp%ntb = 0
   inp%lj1264 = 0
   inp%alpb = 0
   inp%vdwmeth = 0
   inp%ew_type = 0
   inp%inp = 0
   inp%ipb = 0
   inp%igb = 6
   inp%gbsa = 0
   inp%jfastw = 0
   inp%ifqnt = 0
   inp%irism = 0
   inp%rism_verbose = -1
   inp%grdspc1 = 0.8d0
   inp%mdiis_del = 0.7d0
   inp%extdiel = 1.d0
   inp%intdiel = 1.d0
   inp%rgbmax = 25.d0
   inp%saltcon = 0.d0
   inp%cut = 1000.d0
   inp%dielc = 1.d0
   inp%rdt = 0.d0
   inp%ntc = 1
   inp%ntf = 1
   inp%fswitch = -1.d0
   inp%ntr = 0
   inp%ibelly = 0
   inp%mask_from_ref = 0 ! Antoine Marion
   inp%restraint_wt = 0.d0
   inp%restraintmask = ''
   inp%bellymask = ''
   inp%refc = ''

   if (present(gb)) then
      if (gb /= 0 .and. gb /= 1 .and. gb /= 2 .and. gb /= 5 .and. &
          gb /= 6 .and. gb /= 7 .and. gb /= 8 .and. gb /= 10) then
         write(0,*) 'Illegal gb model. Setting to vacuum'
      else
         inp%igb = gb
         inp%extdiel = 78.5
      end if
      ! Vacuum
      if (gb == 0) inp%igb = 6
      ! PB
      if (gb == 10) then
         inp%ipb = 1
         inp%inp = 1
      end if
   end if

end subroutine gas_sander_input

! Sets up input options for a periodic system treated with PME
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
subroutine pme_sander_input(inp)

   use constants, only: NO_INPUT_VALUE

   implicit none

   type(sander_input), intent(out) :: inp

   inp%ntb = 1
   inp%igb = 0
   inp%alpb = 0
   inp%lj1264 = 0
   inp%ipb = 0
   inp%inp = 0
   inp%vdwmeth = 1
   inp%ew_type = 0
   inp%gbsa = 0
   inp%ifqnt = 0
   inp%irism = 0
   inp%rism_verbose = -1
   inp%grdspc1 = 0.8d0
   inp%mdiis_del = 0.7d0
   inp%jfastw = 0
   inp%extdiel = 1.d0
   inp%intdiel = 1.d0
   inp%cut = 8.d0
   inp%rgbmax = 25.d0
   inp%saltcon = 0.d0
   inp%dielc = 1.d0
   inp%rdt = 0.d0
   inp%ntc = 1
   inp%ntf = 1
   inp%fswitch = -1.d0
   inp%ntr = 0
   inp%ibelly = 0
   inp%mask_from_ref = 0 ! Antoine Marion
   inp%restraint_wt = 0.d0
   inp%restraintmask = ''
   inp%bellymask = ''
   inp%refc = ''

end subroutine pme_sander_input

! Initializes a struct with QM options to all default values
!
! Parameters
! ----------
! inp : type(qmmm_input_options)
!     struct of QM input options that will be filled by this subroutine
subroutine qm_sander_input(inp)

   use qmmm_module, only: default_qmmm_input_options
   implicit none

   type(qmmm_input_options), intent(out) :: inp

   call default_qmmm_input_options(inp)

end subroutine qm_sander_input

! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! prmname : string
!     Name of the topology file
! coordinates : double precision(3*natom)
!     Starting coordinates
! inbox : double precision(6)
!     Box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine sander_setup(prmname, coordinates, inbox, input_options, qmmm_options, ierr)

#define USE_PRMTOP_FILE 1
#include "interface_setup.F90"
#undef USE_PRMTOP_FILE

end subroutine sander_setup

#undef rem

subroutine api_mdread1(input_options, ierr)
#define API 1
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#endif /* API */

!$ use omp_lib
   use file_io_dat
   use lmod_driver, only : read_lmod_namelist
   use qmmm_module, only : qmmm_nml, qm_gb
   use constants, only : RETIRED_INPUT_OPTION, zero, one, two, three, seven, &
                         eight, NO_INPUT_VALUE_FLOAT, NO_INPUT_VALUE
   use md_scheme, only: ithermostat, therm_par
   use les_data, only : temp0les
   use stack, only: lastist,lastrst
   use nmr, only: echoin
   use crg_reloc, only: ifcr, cropt, crcut, crskin, crin, crprintcharges
   use sgld, only : isgld, isgsta,isgend,fixcom, &
                    tsgavg,sgft,sgff,sgfd,tempsg,treflf,tsgavp
   use amd_mod, only: iamd,iamdlag,EthreshD,alphaD,EthreshP,alphaP, &
        w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w,igamd
   use scaledMD_mod, only: scaledMD,scaledMD_lambda
   use nbips, only: ips,teips,tvips,teaips,tvaips,raips,mipsx,mipsy,mipsz, &
                    mipso,gridips,dvbips
   use emap,only: temap,gammamap
   use sander_rism_interface, only: rismprm
#ifdef DSSP
   use dssp, only: idssp
#endif /* DSSP */

#ifdef MPI /* SOFT CORE */
   use softcore, only : scalpha,scbeta,ifsc,scmask,logdvdl,dvdl_norest,dynlmb, &
                        sceeorder, tishake
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
   use remd, only  : rem
#endif /* MPI */
   ! Parameter for LIE module
   use linear_response, only: ilrt, lrt_interval, lrtmask
#  ifndef API
   use sander_rism_interface, only: xvvfile, guvfile, huvfile, cuvfile,&
        uuvfile, asympfile, quvFile, chgDistFile, electronMapFile, &
        excessChemicalPotentialfile, solvationEnergyfile, entropyfile, &
        excessChemicalPotentialGFfile, solvationEnergyGFfile, entropyGFfile, &
        excessChemicalPotentialPCPLUSfile, solvationEnergyPCPLUSfile, entropyPCPLUSfile,&
        excessChemicalPotentialUCfile, solvationEnergyUCfile, entropyUCfile,&
        solventPotentialEnergyfile
#  endif /* API */
   use sander_rism_interface, only: rismprm
   use nfe_sander_proxy, only: infe
   implicit none
#  include "box.h"
#  include "def_time.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "../include/md.h"
#  include "../include/memory.h"
#  include "nmr.h"
#  include "tgtmd.h"
#  include "multitmd.h"
#  include "ew_erfc_spline.h"

   character(len=4) watdef(4),watnam,owtnm,hwtnm1,hwtnm2
   _REAL_      dele
   integer     ierr  ! internal error reporting for sander API
   integer     imcdo
   integer     itotst
   integer     inerr
   logical     mdin_cntrl, mdin_lmod, mdin_qmmm  ! true if namelists are in mdin
   integer ::  ifqnt    ! local here --> put into qmmm_nml%ifqnt after read here
   integer     mxgrp
   integer     iemap
   _REAL_      dtemp  ! retired
   _REAL_      dxm  ! retired
   _REAL_      heat  ! retired
   _REAL_      timlim ! retired
   character(11), parameter :: error_hdr =      '| ERROR:   '

! pmemd functionality that is not supported in sander:
   character(len=256) :: scmask1, scmask2, timask1, timask2

#ifdef API
   ! Input options passed in to the setup API routine
   type(sander_input), intent(in) :: input_options
#else
   integer     ifind
   character(len=8) date
   character(len=10) time
   character(len=512) :: char_tmp_512
#endif /* API */

   integer irism, rism_verbose
   double precision grdspc1, mdiis_del
   character(len=8) periodicPotential

!  N.B.: If you make changes to this namelist, you also need to make
!        corresponding changes in ./sander.h and in
!         ../pysander/sander/src/{pysandermodule,pysandermoduletypes}.c

   namelist /cntrl/ irest,ibelly, &
         ntx,ntxo,ntcx,ig,tempi, &
         ntb,ntt,temp0,tautp, &
         ntp,pres0,comp,taup,barostat,mcbarint, &
         nscm,nstlim,t,dt, &
         ntc,ntcc,nconp,tol,ntf,ntn,nsnb, &
         cut,dielc, &
         ntpr,ntwx,ntwv,ntwe,ntwf,ntave,ntpp,ioutfm, &
         ntr,nrc,ntrx,taur,nmropt, &
         ivcap,cutcap,xcap,ycap,zcap,fcap, &
         xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &
         imin,drms,dele,dx0, &
         pencut,ipnlty,iscale,scalm,noeskp, &
         maxcyc,ncyc,ntmin,vlimit, &
         mxsub,ipol,jfastw,watnam,owtnm,hwtnm1,hwtnm2, iesp, &
         skmin, skmax, vv,vfac, tmode, ips, &
         mipsx,mipsy,mipsz,mipso,gridips,raips,dvbips, &
         iamd,iamdlag,EthreshD,alphaD,EthreshP,alphaP, &
         w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w, &
         igamd, &
         ithermostat, therm_par, &
         scaledMD,scaledMD_lambda, &
         iemap,gammamap, &
         isgld,isgsta,isgend,fixcom,tsgavg,sgft,sgff,sgfd,tempsg,treflf,tsgavp,&
         jar, &
         numexchg, repcrd, numwatkeep, hybridgb, reservoir_exchange_step, &
         ntwprt,tausw, &
         ntwr,iyammp,imcdo, &
         plumed,plumedfile, &
         igb,alpb,Arad,rgbmax,saltcon,offset,gbsa,vrand, &
         surften,nrespa,nrespai,gamma_ln,extdiel,intdiel, &
         cut_inner,icfe,clambda,klambda, rbornstat,lastrst,lastist,  &
         itgtmd,tgtrmsd,tgtmdfrc,tgtfitmask,tgtrmsmask, dec_verbose, &
         temp0les,restraintmask,restraint_wt,bellymask, &
         noshakemask,crgmask, &
         mask_from_ref, &
         rdt,ntrelax, &
         ifqnt,ievb, profile_mpi, &
         ipb, inp, &
         gbneckscale, &
         gbalphaH,gbbetaH,gbgammaH, &
         gbalphaC,gbbetaC,gbgammaC, &
         gbalphaN,gbbetaN,gbgammaN, &
         gbalphaOS,gbbetaOS,gbgammaOS, &
         gbalphaP,gbbetaP,gbgammaP, &
         Sh,Sc,Sn,So,Ss,Sp, &
         lj1264, fswitch, &
         ifcr, cropt, crcut, crskin, crin, crprintcharges, &
         ninterface, gamma_ten, infe, baroscalingdir, &
#ifdef MPI /* SOFT CORE */
         scalpha, scbeta, ifsc, scmask, logdvdl, dvdl_norest, dynlmb, &
         sceeorder, &
         ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr, tishake, &
#endif
         ilrt, lrt_interval, lrtmask, &
#ifdef DSSP
         idssp, &
#endif
         irism,&
         vdwmodel, & ! mjhsieh - the model used for van der Waals
         ! retired:
         dtemp, dxm, heat, timlim, &
         ! pmemd functionality that is not supported in sander:
         scmask1, scmask2, timask1, timask2

   ! Define default water residue name and the names of water oxygen & hydrogens

   data watdef/'WAT ','O   ','H1  ','H2  '/

   !     ----- READ THE CONTROL DATA AND OPEN DIFFERENT FILES -----

#ifndef API /* If this is NOT the API */
   if (mdout /= "stdout" ) &
         call amopen(6,mdout,owrite,'F','W')
   call amopen(5,mdin,'O','F','R')
   write(6,9308)
   call date_and_time( DATE=date, TIME=time )
   write(6,'(12(a))') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)

   ! Write the path of the current executable and working directory
   call get_command_argument(0, char_tmp_512)
   write(6,'(/,a,a)') '|   Executable path: ', trim(char_tmp_512)
   call getcwd(char_tmp_512)
   write(6,'(a,a)') '| Working directory: ', trim(char_tmp_512)
! Write the hostname if we can get it from environment variable
! Note: get_environment_variable is part of the F2003 standard but seems
!       to be supported by GNU, Intel, IBM and Portland (2010+) compilers
   call get_environment_variable("HOSTNAME", char_tmp_512, inerr)
   if (inerr .eq. 0) then
     write(6,'(a,a,/)') '|          Hostname: Unknown'
   else
     write(6,'(a,a,/)') '|          Hostname: ', trim(char_tmp_512)
   end if

   if (owrite /= 'N') write(6, '(2x,a)') '[-O]verwriting output'

   ! Echo the file assignments to the user:

   write(6,9700) 'MDIN'   ,mdin(1:70)  , 'MDOUT' ,mdout(1:70) , &
         'INPCRD' ,inpcrd(1:70), 'PARM'  ,parm(1:70)  , &
         'RESTRT',restrt(1:70) , 'REFC'  ,refc(1:70)  , &
         'MDVEL' ,mdvel(1:70)  , 'MDFRC' ,mdfrc(1:70) , &
         'MDEN'   ,mden(1:70)  , &
         'MDCRD' ,mdcrd(1:70)  , 'MDINFO' ,mdinfo(1:70), &
         'MTMD'  ,mtmd(1:70)   , 'INPDIP', inpdip(1:70), &
         'RSTDIP', rstdip(1:70), 'INPTRAJ', inptraj(1:70)
#  ifdef MPI
   write(6,9702) 'REMLOG',     trim(remlog), &
                 'REMTYPE',    trim(remtype), &
                 'REMSTRIP',   trim(remstripcoord), &
                 'SAVEENE',    trim(saveenefile), &
                 'CLUSTERINF', trim(clusterinfofile), &
                 'RESERVOIR',  trim(reservoirname), &
                 'REMDDIM',    trim(remd_dimension_file)
#  endif
   if (len_trim(xvvfile) > 0) &
        write(6,9701) 'Xvv', trim(xvvfile)
   if (len_trim(guvfile) > 0) &
        write(6,9701) 'Guv', trim(Guvfile)
   if (len_trim(huvfile) > 0) &
        write(6,9701) 'Huv', trim(Huvfile)
   if (len_trim(cuvfile) > 0) &
        write(6,9701) 'Cuv', trim(Cuvfile)
   if (len_trim(uuvfile) > 0) &
        write(6,9701) 'Uuv', trim(Uuvfile)
   if (len_trim(asympfile) > 0) &
        write(6,9701) 'Asymptotics', trim(asympfile)
   if (len_trim(quvfile) > 0) &
        write(6,9701) 'Quv', trim(Quvfile)
   if (len_trim(chgDistfile) > 0) &
        write(6,9701) 'ChgDist', trim(chgDistfile)
   if (len_trim(electronMapFile) > 0) &
        write(6,9701) 'ElectronMap', trim(electronMapFile)
   if (len_trim(excessChemicalPotentialfile) > 0) &
        write(6,9701) 'ExChem', trim(excessChemicalPotentialfile)
   if (len_trim(solvationEnergyfile) > 0) &
        write(6,9701) 'SolvEne', trim(solvationEnergyfile)
   if (len_trim(entropyfile) > 0) &
        write(6,9701) 'Entropy', trim(entropyfile)
   if (len_trim(excessChemicalPotentialGFfile) > 0) &
        write(6,9701) 'ExChGF', trim(excessChemicalPotentialGFfile)
   if (len_trim(solvationEnergyGFfile) > 0) &
        write(6,9701) 'SolvEneGF', trim(solvationEnergyGFfile)
   if (len_trim(entropyGFfile) > 0) &
        write(6,9701) '-TS_GF', trim(entropyGFfile)
   if (len_trim(excessChemicalPotentialPCPLUSfile) > 0) &
        write(6,9701) 'exchemPCPLUS', trim(excessChemicalPotentialPCPLUSfile)
   if (len_trim(solvationEnergyPCPLUSfile) > 0) &
        write(6,9701) 'SolvEnePCPLUS', trim(solvationEnergyPCPLUSfile)
   if (len_trim(entropyPCPLUSfile) > 0) &
        write(6,9701) '-TS_PCPLUS', trim(entropyPCPLUSfile)
   if (len_trim(excessChemicalPotentialUCfile) > 0) &
        write(6,9701) 'exchemUC', trim(excessChemicalPotentialUCfile)
   if (len_trim(solvationEnergyUCfile) > 0) &
        write(6,9701) 'SolvEneUC', trim(solvationEnergyUCfile)
   if (len_trim(entropyUCfile) > 0) &
        write(6,9701) '-TS_UC', trim(entropyUCfile)
   if (len_trim(solventPotentialEnergyfile) > 0) &
        write(6,9701) 'PotUV', trim(solventPotentialEnergyfile)

   ! Echo the input file to the user:
   call echoin(5,6)
   !     ----- READ DATA CHARACTERIZING THE MD-RUN -----
   read(5,'(a80)') title
   !       ----read input in namelist format, first setting up defaults
#endif /* ifndef API */

   ierr = 0

   dtemp = RETIRED_INPUT_OPTION
   dxm   = RETIRED_INPUT_OPTION
   heat  = RETIRED_INPUT_OPTION
   timlim = RETIRED_INPUT_OPTION
! pmemd functionality that is not supported in sander:
   scmask1 = ''
   scmask2 = ''
   timask1 = ''
   timask2 = ''

   irest = 0
   ibelly = 0
   mask_from_ref = 0 ! Antoine Marion
   ipol = RETIRED_INPUT_OPTION
   iesp = 0
   ntx = 1
   ntxo = NO_INPUT_VALUE
   ig = 71277
   tempi = ZERO
   ntb = NO_INPUT_VALUE
   ntt = 0
   temp0 = 300.0d0
! MIDDLE SCHEME{ 
   ithermostat = 1
   therm_par = 5.0d0
! } 
! PLUMED
   plumed = 0
   plumedfile = 'plumed.dat'
! END PLUMED
#ifdef LES
   ! alternate temp for LES copies, if negative then use single bath
   ! single bath not the same as 2 baths with same target T
   temp0les = -ONE
#endif
   rdt = 0

   tautp = ONE
   ntp = 0
   fswitch = -1.d0
   barostat = 1
   mcbarint = 100
   pres0 = ONE
   comp = 44.6d0
   taup = ONE
   npscal = 1
   nscm = 1000
   nstlim = 1
   t = ZERO
   dt = 0.001d0
   ntc = 1
   tol = 0.00001
   ntf = 1
   nsnb = 25
   cut =  NO_INPUT_VALUE_FLOAT
   dielc = ONE
   ntpr = 50
!RCW - Amber 16 change ntwr to equal nstlim
!0 here is to detect if user set it or not.
   ntwr = 0
   ntwx = 0
   ntwv = 0
   ntwf = 0
   ntwe = 0
   ipb = 0
   inp = 2
   irism = 0
   rism_verbose = -1
   grdspc1 = 0.8d0
   mdiis_del = 0.7d0
   ntave = 0
!RCW: Amber 16 default to netcdf if support is compiled in.
   ioutfm = 1
   ntr = 0
   ntrx = 1
   ivcap = 0
   natcap = 0
   fcap = 1.5d0
   cutcap = 0.0d0
   xcap = 0.0d0
   ycap = 0.0d0
   zcap = 0.0d0
   forth = 1.5d0
   xlorth = -1.0d0
   ylorth = -1.0d0
   zlorth = -1.0d0
   xorth = 47114711.0d0
   yorth = 47114711.0d0
   zorth = 47114711.0d0
   numexchg = 0
   repcrd   = 1
   lj1264 = -1

   profile_mpi = 0 !whether to write profile_mpi timing file - default = 0 (NO).

   ! number of waters to keep for hybrid model,
   ! numwatkeep: the number of closest
   ! waters to keep. close is defined as close to non-water.
   ! for simulations with ions, ions should be stripped too
   ! or at least ignored in the "closest" calculation. this
   ! is not currently done.

   ! if it stays at -1 then we keep all waters
   ! 0 would mean to strip them all

    numwatkeep=-1

   ! hybridgb: gb model to use with hybrid REMD.
   hybridgb=0

   ! Reservoir Exchange Frequency
   reservoir_exchange_step=2

   ! carlos targeted MD, like ntr

   itgtmd=0
   tgtrmsd=0.
   tgtmdfrc=0.
   tgtfitmask=''
   tgtrmsmask=''

   pencut = 0.1d0
   taumet = 0.0001d0
   omega = 500.0d0
   ipnlty = 1
   scalm = 100.0d0
   iscale = 0
   noeskp = 1
   nmropt = 0
   jar = 0
   tausw = 0.1d0
   imin = 0
   isftrp = 0
   rwell = ONE
   maxcyc = 1
   ncyc = 10
   ntmin = 3
   dx0 = 0.01d0
   drms = 1.0d-4
   vlimit = 20.0d0
   mxsub = 1
   jfastw = 0
   watnam = '    '
   owtnm =  '    '
   hwtnm1 = '    '
   hwtnm2 = '    '
   ntwprt = 0
   igb = 0
   alpb = 0
   Arad = 15.0d0
   rgbmax = 25.d0
   saltcon = ZERO

   !  default offset depends on igb value, and users need to
   !  be able to modify it, so we need to set a dummy value. if it's still the
   !  dummy after we read the namelist, we set the default based on igb. if not,
   !  we leave it at what the user set.
   !  best solution would be to create a GB namelist.
   offset = -999999.d0
   gbneckscale = -999999.d0

   iyammp = 0
   imcdo = -1
   gbsa = 0
   vrand=1000
   surften = 0.005d0
   nrespa = 1
   nrespai = 1
   irespa = 1
   gamma_ln = ZERO
   extdiel = 78.5d0
   intdiel = ONE
   gbgamma = ZERO
   gbbeta = ZERO
   gbalpha = ONE

   ! set default parameters for igb = 8
   ! NOTE THAT NONE OF THESE ARE USED UNLESS IGB=8, SO USERS SHOULD NOT EVEN SET
   ! THEM
   gbalphaH = 0.788440d0
   gbbetaH = 0.798699d0
   gbgammaH = 0.437334d0
   gbalphaC = 0.733756d0
   gbbetaC = 0.506378d0
   gbgammaC = 0.205844d0
   gbalphaN = 0.503364d0
   gbbetaN = 0.316828d0
   gbgammaN = 0.192915d0
   gbalphaOS = 0.867814d0
   gbbetaOS = 0.876635d0
   gbgammaOS = 0.387882d0
   gbalphaP = 1.0d0    !P parameters are not optimized yet
   gbbetaP = 0.8d0     !P parameters are not optimized yet
   gbgammaP = 4.85d0   !P parameters are not optimized yet
   !scaling parameters below will only be used for igb=8.
   ! the actual code does not use these variables, it uses X(l96)
   ! if igb=8, we will use these to set the X(l96) array.
   Sh = 1.425952d0
   Sc = 1.058554d0
   Sn = 0.733599d0
   So = 1.061039d0
   Ss = -0.703469d0
   Sp = 0.5d0          !P parameters are not optimized for protein
   ! update gbneck2nu pars
   ! using offset and gb_neckscale parameters from GB8-protein for gbneck2nu
   ! Scaling factors
   ! name of variables are different from sander's igb8
   ! we use below names in pmemd
   ! gbneck2nu
   screen_hnu = 1.696538d0
   screen_cnu = 1.268902d0
   screen_nnu = 1.4259728d0
   screen_onu = 0.1840098d0
   screen_pnu = 1.5450597d0
   !alpha, beta, gamma for each atome element
   gb_alpha_hnu = 0.537050d0
   gb_beta_hnu = 0.362861d0
   gb_gamma_hnu = 0.116704d0
   gb_alpha_cnu = 0.331670d0
   gb_beta_cnu = 0.196842d0
   gb_gamma_cnu = 0.093422d0
   gb_alpha_nnu = 0.686311d0
   gb_beta_nnu = 0.463189d0
   gb_gamma_nnu = 0.138722d0
   gb_alpha_osnu = 0.606344d0
   gb_beta_osnu = 0.463006d0
   gb_gamma_osnu = 0.142262d0
   gb_alpha_pnu = 0.418365d0
   gb_beta_pnu = 0.290054d0
   gb_gamma_pnu = 0.1064245d0
   ! End gbneck2nu

   iconstreff = 0
   cut_inner = EIGHT
   icfe = 0
   clambda = ZERO
   klambda = 1
   ievb = 0
   rbornstat = 0
   ! added a flag to control output of BDC/SDC synonymous with MMPBSA.py's
   ! version of the same variable.
   dec_verbose = 3
   lastrst = 1
   lastist = 1
   restraintmask=''
   restraint_wt = ZERO
   bellymask=''
   noshakemask=''
   crgmask=''

   ntrelax = 500 ! how long to let waters relax
   skmin = 50 !used by neb calculation
   skmax = 100 !used by neb calculation
   vv = 0 !velocity verlet -- off if vv/=1
   vfac = 0 !velocity verlet scaling factor, 0 by default
   tmode = 1 !default tangent mode for NEB calculation

   ifqnt = NO_INPUT_VALUE
   rism_verbose = NO_INPUT_VALUE
   irism = NO_INPUT_VALUE
   grdspc1 = NO_INPUT_VALUE_FLOAT
   mdiis_del = NO_INPUT_VALUE_FLOAT

   ifcr = 0 ! no charge relocation
   cropt = 0 ! 1-4 EEL is calculated with the original charges
   crcut = 3.0
   crskin = 2.0
   crin = ''
   crprintcharges = 0

   ips = 0    ! no isotropic periodic sum
   raips=-1.0d0   ! automatically determined
   mipsx=-1   ! number of grids in x direction, <0 for automatically determined
   mipsy=-1   ! number of grids in y direction, <0 for automatically determined
   mipsz=-1   ! number of grids in z direction, <0 for automatically determined
   mipso=4    ! default 4th order b-spline
   gridips=2   ! grid size. used to determine grid number if not defined
   dvbips=1.0d-8   ! Volume change tolerance. aips will be done when change more than dvbips

   iamd = 0 ! No accelerated MD used
   iamdlag = 0 !frequency of boosting in steps
   EthreshD = 0.d0
   alphaD = 0.d0
   EthreshP = 0.d0
   alphaP = 0.d0
   w_amd = 0 ! windowed amd
   EthreshD_w = 0.d0
   alphaD_w = 0.d0
   EthreshP_w = 0.d0
   alphaP_w = 0.d0

   igamd = 0 ! No GaMD used

   scaledMD = 0 ! No scaled MD used
   scaledMD_lambda = 0.d0

   iemap=0     ! no emap constraint
   gammamap=1     ! default friction constant for map motion, 1/ps
   isgld = 0   ! no self-guiding
   isgsta=1    ! Begining index of SGLD range
   isgend=0    ! Ending index of SGLD range
   fixcom=-1    ! fix center of mass in SGLD simulation
   tsgavg=0.2d0    !  Local averaging time of SGLD simulation
   sgft=-1.0d3      !  Guiding factor of SGLD simulation
   sgff=-1.0d3      !  Guiding factor of SGLD simulation
   sgfd=-1.0d3      !  Guiding factor of SGLD simulation
   tempsg=0.0d0    !  Guiding temperature of SGLD simulation
   treflf=0.0d0    !  Reference low frequency temperature of SGLD simulation
   tsgavp=2.0d0    !  Convergency time of SGLD simulation

   !     Check to see if "cntrl" namelist has been defined.
   mdin_cntrl=.false.
   mdin_qmmm = .false.
   mdin_ewald=.false.
   mdin_pb=.false.
   mdin_lmod=.false.
#ifdef MPI /* SOFT CORE */
   scalpha=0.5
   scbeta=12.0
   sceeorder=2
   ifsc=0
   logdvdl=0
   dvdl_norest=0
   dynlmb=0.0
   ifmbar=0
   bar_intervall=100
   bar_l_min=0.1
   bar_l_max=0.9
   bar_l_incr=0.1
   tishake = 0
#endif
   ilrt = 0
   lrt_interval = 50
   lrtmask=''
#ifdef DSSP
   idssp = 0
#endif

!  ntp = 2 anisotropic directional pressure scaling
   baroscalingdir = 0     !default=0, random direction pressure scaling per step

#ifdef API
   igb = input_options%igb
   alpb = input_options%alpb
   gbsa = input_options%gbsa
   lj1264 = input_options%lj1264
   ipb = input_options%ipb
   inp = input_options%inp
   vdwmeth = input_options%vdwmeth
   ew_type = input_options%ew_type
   extdiel = input_options%extdiel
   intdiel = input_options%intdiel
   rgbmax = input_options%rgbmax
   saltcon = input_options%saltcon
   cut = input_options%cut
   dielc = input_options%dielc
   ifqnt = input_options%ifqnt
   irism = input_options%irism
   rism_verbose = input_options%rism_verbose
   grdspc1 = input_options%grdspc1
   mdiis_del = input_options%mdiis_del
   jfastw = input_options%jfastw
   ntf = input_options%ntf
   ntc = input_options%ntc
   fswitch = input_options%fswitch
   ntr = input_options%ntr
   ibelly = input_options%ibelly
   restraint_wt = input_options%restraint_wt
   restraintmask = input_options%restraintmask
   bellymask = input_options%bellymask
   refc = input_options%refc
   mask_from_ref = input_options%mask_from_ref ! Antoine Marion
#  ifdef LES
   rdt = input_options%rdt
#  endif
#else /* NOT the API */
   call nmlsrc('cntrl',5,ifind)
   if (ifind /= 0) mdin_cntrl=.true.
   rewind(5)

   call nmlsrc('ewald',5,ifind)
   if (ifind /= 0) mdin_ewald=.true.
   rewind 5

   call nmlsrc('pb',5,ifind)
   if (ifind /= 0) mdin_pb=.true.
   rewind 5

   call nmlsrc('qmmm', 5, ifind)
   if (ifind /= 0) mdin_qmmm = .true.
   rewind 5

   call nmlsrc('lmod',5,ifind)
   if (ifind /= 0) mdin_lmod=.true.
   rewind 5

   rewind 5
   if ( mdin_cntrl ) then
      read(5,nml=cntrl,err=999)
   else
      write(6, '(1x,a,/)') 'Could not find cntrl namelist'
      FATAL_ERROR
   end if
#endif /* ifdef API */

   ! Set ntwr to nstlim if the user did not specify anything.
   ! Note nstlim can actually be set to zero but zero is not permitted for ntwr
   ! so if it is still zero set it to 1.
   if ( nstlim .eq. 0 ) then
      ntwr = 1
   else
      if ( ntwr .eq. 0 ) then
         if ( numexchg .eq. 0 ) then
            ntwr = nstlim
         else
            ntwr = nstlim * numexchg   ! for replica exchange
         end if
      end if
   end if

   if ( ntwr .eq. 0 ) ntwr = nstlim
   if ( ntwr .eq. 0 ) ntwr = 1

   if ( igb == 10 .and. ipb == 0 ) ipb = 2
   if ( igb == 0  .and. ipb /= 0 ) igb = 10

   if (plumed.eq.1) then
     write(6, '(1x,a,/)') 'PLUMED is on'
     write(6, '(1x,a,a,/)') 'PLUMEDfile is ',plumedfile
   endif

   if (ifqnt == NO_INPUT_VALUE) then
      ifqnt = 0 ! default value
      if (mdin_qmmm) then
         write(6, '(1x,a,/)') &
            '| WARNING qmmm namelist found, but ifqnt was not set! QMMM NOT &
            &active.'
      end if
   end if

   if (irism == NO_INPUT_VALUE) then
      irism = 0 ! default value
   end if
   if (rism_verbose == NO_INPUT_VALUE) then
      rism_verbose = -1 ! default value
   end if
   if (grdspc1 == NO_INPUT_VALUE) then
      grdspc1 = 0.8d0 ! default value
   end if
   if (mdiis_del == NO_INPUT_VALUE) then
      mdiis_del = 0.7d0 ! default value
   end if

   ! middle scheme is requested {
      if (ithermostat < 0 .or. ithermostat > 2) then
         write(6,'(1x,a,/)') &
            'Middle scheme: ithermostat is only available for 0-2'
         FATAL_ERROR
      end if
      if (therm_par < 0d0) then
         write(6,'(1x,a,/)') 'Middle scheme: therm_par MUST be non-negative'
         FATAL_ERROR
      end if
   ! }

   ! Now that we've read the input file, set up the defaults for variables
   ! whose values depend on other input values (ntb, cut)
   if (ntb == NO_INPUT_VALUE) then
      if (ntp > 0) then
         ntb = 2
      else if (igb > 0) then
         ntb = 0
      else
         ntb = 1
      end if
   end if

   if (ntxo == NO_INPUT_VALUE) then
      ntxo = 2
   end if

   if (cut == NO_INPUT_VALUE_FLOAT) then
      if (igb == 0) then
         cut = EIGHT
      else
         cut = 9999.d0
      end if
   end if

   ! Force igb=6 to get vacuum electrostatics or igb=0 for periodic
   ! boundary conditions. This must be done ASAP to ensure SANDER's
   ! electrostatics are initialized properly.

   rismprm%rism=irism
   rismprm%verbose=rism_verbose
   rismprm%grdspc(:)=grdspc1
   rismprm%mdiis_del=mdiis_del

   if (irism /= 0) then
      periodicPotential = 'pme'
      igb = 0
   end if

   if (ifqnt>0) then
      qmmm_nml%ifqnt = .true.
      if (saltcon /= 0.0d0) then
         qm_gb%saltcon_on = .true.
      else
         qm_gb%saltcon_on = .false.
      end if
      if (alpb == 1) then
         qm_gb%alpb_on = .true.
      else
         qm_gb%alpb_on = .false.
      end if
      if (igb == 10 .or. ipb /= 0) then
         write(6, '(1x,a,/)') 'QMMM is not compatible with Poisson Boltzmann (igb=10 or ipb/=0).'
         FATAL_ERROR
      end if
   else
      qmmm_nml%ifqnt = .false.
   end if

   if ( mdin_lmod ) then
      rewind 5
      call read_lmod_namelist()
   end if

   !--------------------------------------------------------------------
   !     --- vars have been read ---
   !--------------------------------------------------------------------

#ifndef API
   write(6,9309)
#endif

   ! emit warnings for retired cntrl namelist variables

   if ( dtemp /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: dtemp has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( dxm /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: dxm has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
            ! '  The step length will be unlimited.'
   end if
   if ( heat /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: heat has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if
   if ( timlim /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/)') 'Warning: timlim has been retired.', &
            '  Check the Retired Namelist Variables Appendix in the manual.'
   end if

   ! emit warnings for pmemd functionality that is not supported in sander:

   if ( scmask1 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: scmask1 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Softcore Potentials section in the manual.'
   end if
   if ( scmask2 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: scmask2 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Softcore Potentials section in the manual.'
   end if
   if ( timask1 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: timask1 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Thermodynamic Integration section in the manual.'
   end if
   if ( timask2 /= '' ) then
      write(6,'(/,a,/,a,/,a/)') 'Warning: timask2 is a pmemd variable.', &
            '  It is not supported in sander.', &
            '  Check the Thermodynamic Integration section in the manual.'
   end if

! baroscalingdir valid options
   if (baroscalingdir > 0) then
      if (baroscalingdir < 0 .or. baroscalingdir > 3) then
         write(6,'(/2x,a)') &
         'Invalid baroscalingdir value. baroscalingdir must be between 0 and 3'
         FATAL_ERROR
      end if
      if (ntb /= 2) then
         write(6,'(/2x,a)') &
         'ntb invalid. ntb must be 2 for directional pressure scaling.'
         FATAL_ERROR
      end if
      if (ntp /= 2) then
         write(6,'(/2x,a)') &
         'ntp invalid. ntp must be 2 for directional pressure scaling.'
         FATAL_ERROR
      end if
   end if

! MC Barostat valid options. Some of these may work, but disable them until they
! are fully tested.

   if (ntp > 0 .and. barostat == 2) then
      inerr = 0
      if (icfe /= 0) then
         write(6, '(/2x,a)') 'TI is not compatible with the MC Barostat'
         inerr = 1
      end if
#ifdef LES
      write(6, '(/2x,a)') 'LES is not compatible with the MC Barostat'
      inerr = 1
#endif
      ! Any others? Hopefully most or all of the above can be made compatible.
      if (inerr == 1) then
         FATAL_ERROR
      end if
   end if

#ifndef API
   call printflags()
#endif

   !--------------------------------------------------------------------
   ! If user has requested ewald electrostatics, read some more input
   !--------------------------------------------------------------------

#ifdef API
   if( igb == 0 .and. ipb == 0 ) call load_ewald_info(ntp)
#else
   if( igb == 0 .and. ipb == 0 ) call load_ewald_info(inpcrd,ntp)
#endif

   !--------------------------------------------------------------------
   ! parameters for IPS and for SGLD:
   ! ips=1  3D IPS for electrostatic and Lennard-Jones potentials
   ! ips=2  3D IPS for electrostatic potential only
   ! ips=3  3D IPS for  Lennard-Jones potential only
   ! ips=4  3D IPS/DFFT for electrostatic and Lennard-Jones potentials
   ! ips=5  3D IPS/DFFT for electrostatic potential only
   ! ips=6  3D IPS/DFFT for Lennard-Jones potential only
   !--------------------------------------------------------------------

   teips=.false.
   tvips=.false.
   teaips=.false.
   tvaips=.false.
   if((ips-4)*(ips-6) == 0 )tvaips =.true.
   if ( (ips-4)*(ips-5) == 0 )teaips =.true.
   if( tvaips.OR.( (ips -1)*(ips-3) == 0 ))tvips =.true.
   if( teaips.OR.((ips -1)*(ips-2) == 0 ))teips =.true.
   if( teips ) then
      use_pme = 0
      eedmeth = 6
   end if
   if( tvips ) then
      vdwmeth = 2
      if(use_pme/=0.and.tvaips)then
        mipsx=nfft1   ! number of grids in x direction, <0 for automatically determined
        mipsy=nfft2   ! number of grids in y direction, <0 for automatically determined
        mipsz=nfft3   ! number of grids in z direction, <0 for automatically determined
        mipso=order    ! default 6th order b-spline
      endif
   end if
   temap=iemap>0
   ishake = 0
   if (ntc > 1) ishake = 1

   !--------------------------------------------------------------------
   ! Set up some parameters for AMD simulations:
   ! AMD initialization
   ! iamd=0 no boost is used, 1 boost on the total energy,
   ! 2 boost on the dohedrals, 3 boost on dihedrals and total energy
   !--------------------------------------------------------------------
   if(iamd.gt.0)then
      if(iamd.eq.1)then !only total potential energy will be boosted
         EthreshD=0.d0
         alphaD=0.d0
      else if(iamd.eq.2)then !only dihedral energy will be boosted
         EthreshP=0.d0
         alphaP=0.d0
      endif
      if(w_amd.gt.0)then
         if(iamd.eq.1)then !only total potential energy will be boosted
            EthreshD_w=0.d0
            alphaD_w=0.d0
         else if(iamd.eq.2)then !only dihedral energy will be boosted
            EthreshP_w=0.d0
            alphaP_w=0.d0
         endif
#ifndef API
         write(6,'(a,i3)')'| Using Windowed Accelerated MD (wAMD) &
                          &LOWERING BARRIERS to enhance sampling w_amd =', w_amd
         write(6,'(a,2f22.12)')'| AMD boost to total energy: EthreshP,alphaP',&
                           EthreshP, alphaP
         write(6,'(a,2f22.12)')'| AMD boost to dihedrals: EthreshD,alphaD',&
                           EthreshD,alphaD
         write(6,'(a,2f22.12)')'| AMD extra parameters boost to total energy: &
                           &EthreshP_w,alphaP_w', EthreshP_w, alphaP_w
         write(6,'(a,2f22.12)')'| AMD extra parameters boost to dihedrals: &
                           &EthreshD_w,alphaD_w', EthreshD_w, alphaD_w
      else
         write(6,'(a,i3)')'| Using Accelerated MD (AMD) RASING VALLEYS to &
                           &enhance sampling iamd =',iamd
         write(6,'(a,2f22.12)')'| AMD boost to total energy: EthreshP,alphaP', &
                           EthreshP, alphaP
         write(6,'(a,2f22.12)')'| AMD boost to dihedrals: EthreshD,alphaD', &
                           EthreshD, alphaD
#endif
      endif
   endif


   !--------------------------------------------------------------------
   ! Set up some parameters for scaledMD simulations:
   ! scaledMD initialization
   ! scaledMD=0 no scaling is used, 1 scale the potential energy
   !--------------------------------------------------------------------
#ifndef API
   if(scaledMD.gt.0)then
      write(6,'(a,i3)')'| Using Scaled MD to enhance sampling scaledMD =',&
                        scaledMD
      write(6,'(a,f22.12)')'| scaledMD scaling factor lambda: ',scaledMD_lambda
   endif
#endif

   !--------------------------------------------------------------------
   ! Set up some parameters for GB simulations:
   !--------------------------------------------------------------------
   ! update offset = 0.09d0 for igb /= 8
   ! I add this step because I want to use different offset value as default value
   ! for igb = 8

   if ( igb == 8 ) then
     if (offset == -999999.d0) then
        offset = 0.195141d0  !set to default for igb=8
     end if
     if (gbneckscale == -999999.d0) then
        gbneckscale = 0.826836d0
     end if
   else
      ! not igb=8, use old defaults
      if (offset == -999999.d0) then
         offset = 0.09d0
      end if
      if (gbneckscale == -999999.d0) then
         gbneckscale = 0.361825d0
      end if
   endif

   if( igb == 2 .or. hybridgb == 2 ) then
      !       --- use our best guesses for Onufriev/Case GB  (GB^OBC I)

      gbgamma = 2.90912499999d0  ! (the "99999" to force roundoff on print)
      gbbeta = ZERO
      gbalpha = 0.8d0
   end if

   if( igb == 5 .or. hybridgb == 5 ) then

      !       --- use our second best guesses for Onufriev/Case GB (GB^OBC II)

      gbgamma = 4.850d0
      gbbeta = 0.8d0
      gbalpha = ONE
   end if

   if( igb == 7 ) then

      !       --- use parameters for Mongan et al. CFA GBNECK

      gbgamma = 2.50798245d0
      gbbeta = 1.90792938d0
      gbalpha = 1.09511284d0
   end if

   ! -------------------------------------------------------------------
   ! If the user has requested NMR restraints, do a cursory read of the
   ! restraints file(s) now to determine the amount of memory necessary
   ! for these restraints:
   ! -------------------------------------------------------------------

   if (jar == 1 ) nmropt = 1
   intreq = 0
   irlreq = 0
   if (nmropt > 0) then
      mxgrp = 0
      itotst = 1

      ! Set ITOTST to 0 if IMIN equals 1 (i.e. if minimization, not dynamics)
      ! This will cause any "time-averaged" requests to be over-ridden.

      if (imin == 1) then
         itotst = 0
      end if
      !         CALL AMOPEN(31,NMR,'O','F','R')
      call restlx(5,itotst,mxgrp,dt,6,ierr)
      !         CLOSE(31)
   end if

   ! Set the definition of the water molecule. The default definition is in
   ! WATDEF(4).

   read(watdef(1),'(A4)') iwtnm
   read(watdef(2),'(A4)') iowtnm
   read(watdef(3),'(A4)') ihwtnm(1)
   read(watdef(4),'(A4)') ihwtnm(2)
   if (watnam /= '    ') read(watnam,'(A4)') iwtnm
   if (owtnm /= '    ') read(owtnm, '(A4)') iowtnm
   if (hwtnm1 /= '    ') read(hwtnm1,'(A4)') ihwtnm(1)
   if (hwtnm2 /= '    ') read(hwtnm2,'(A4)') ihwtnm(2)

   return

#ifndef API
   999 continue   ! bad cntrl read
   write(6,*) 'error in reading namelist cntrl'
   FATAL_ERROR

   ! --- input file polar opts read err trapping:

   9308 format(/10x,55('-'),/10x, &
         'Amber 18 SANDER                              2018', &
         /10x,55('-')/)
   9309 format(/80('-')/'   1.  RESOURCE   USE: ',/80('-')/)
   9700 format(/,'File Assignments:',/,15('|',a6,': ',a,/))
   9701 format('|',a6,': ',a)
#  ifdef MPI
   9702 format(7('|',a10,': ',a,/))
#  endif /* MPI */
#endif /* API */

#undef FATAL_ERROR
#undef API
end subroutine api_mdread1

subroutine api_mdread2(x, ix, ih, ierr)
#define API 1
#ifdef API
#  define FATAL_ERROR ierr = 1; return
#  define DELAYED_ERROR ierr = 1; return
#else
#  define FATAL_ERROR call mexit(6, 1)
#  define DELAYED_ERROR inerr = 1
#endif /* API */

   use lmod_driver, only : LMOD_NTMIN_LMOD, LMOD_NTMIN_XMIN, write_lmod_namelist
   use findmask
#ifdef LES
   use genbornles, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#else
   use genborn, only: isnucat ! gbneck2nu: check if atom belongs to nuc or protein
#endif
   use constants, only : ZERO, ONE, TWO
   use parms, only: req
   use nbips, only: ips
   use amd_mod, only: iamd,EthreshD,alphaD,EthreshP,alphaP, &
        w_amd,EthreshD_w,alphaD_w,EthreshP_w,alphaP_w,igamd
   use nblist, only: a,b,c,alpha,beta,gamma,nbflag,skinnb,sphere,nbtell,cutoffnb
   use md_scheme, only: therm_par
   use file_io_dat
   use sander_lib, only: upper
#ifdef LES
   use les_data, only : lestyp, lesfac, nlesty, temp0les
#endif
#ifndef API
   use emap, only : temap, emap_options
   use qmmm_module, only : qmmm_nml, qmmm_vsolv, qmmm_struct
   use qmmm_vsolv_module, only : print
   use linear_response, only : lrt_interval
   use crg_reloc, only: ifcr, cropt, crcut, crskin, crprintcharges
#else
#  ifdef LES
   use qmmm_module, only : qmmm_nml
#  endif
#endif /* API */
   use qmmm_module, only: get_atomic_number, qm_gb
#ifdef MPI
   use softcore, only : ifsc, scalpha, scbeta, dvdl_norest, &
                        sceeorder, logdvdl, dynlmb
   use mbar, only : ifmbar, bar_intervall, bar_l_min, bar_l_max, bar_l_incr
! REMD
   use remd, only : rem, rremd
   use sgld, only : isgld ! for RXSGLD
#endif /* MPI */
   use linear_response, only: ilrt, lrtmask
   use nfe_sander_proxy, only: infe
   use commandline_module, only: cpein_specified

   implicit none
   _REAL_ x(*)
#  include "../include/memory.h"
   integer ix(lasti)
   character(len=4) ih(*)
   integer ierr, nbond
   integer atom1,atom2
   logical belly,konst
   character(len=2) atype
   integer atomicnumber, hybridization
   integer ngrp,inerr,nr,ir,i,mxresat,j
   integer noshakegp( natom ), natnos
   integer ier
   logical errFlag
   _REAL_ emtmd, wallc
#ifndef LES
   logical newstyle
#endif /* LES */
#ifndef API
   integer ntmp
#endif /* API */

#ifdef MPI
   !     =========================== AMBER/MPI ===========================
#  ifdef MPI_DOUBLE_PRECISION
#     undef MPI_DOUBLE_PRECISION
#  endif
   include 'mpif.h'
#  include "parallel.h"
   integer ist(MPI_STATUS_SIZE), partner, nbonh_c, num_noshake_c
   integer nquant_c, noshake_overlap_c
   integer crggp( natom )
   _REAL_ val
   _REAL_, allocatable :: repvals(:)
   !     ========================= END AMBER/MPI =========================
#endif /* MPI */
#  include "../include/md.h"
#  include "box.h"
#  include "nmr.h"
#  include "extra_pts.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "ew_mpole.h"
#  include "ew_legal.h"
#  include "def_time.h"
#  include "tgtmd.h"
#  include "multitmd.h"
   character(11), parameter :: error_hdr =      '| ERROR:   '

   ! -------------------------------------------------------------------
   !      --- set up resat array, containing string identifying
   !          residue for each atom
   ! -------------------------------------------------------------------
   ierr = 0
   mxresat = min( natom, matom )
   ir = 0
   do i=1,mxresat
      if (i >= ix(ir+i02)) ir=ir+1
      write(resat(i),'(a4,1x,a4,i4)') ih(m04+i-1), &
            ih(ir+m02-1),ir
      !                                    ---null terminator:
      resat(i)(14:14) = char(0)
   end do
   close(unit=8)

   ! -------------------------------------------------------------------
   !     ----- SET THE DEFAULT VALUES FOR SOME VARIABLES -----
   ! -------------------------------------------------------------------

   nrp = natom

#ifndef API
   if (ifbox == 1) write(6, '(/5x,''BOX TYPE: RECTILINEAR'')')
   if (ifbox == 2) write(6, '(/5x,''BOX TYPE: TRUNCATED OCTAHEDRON'')')
   if (ifbox == 3) write(6, '(/5x,''BOX TYPE: GENERAL'')')
#endif

   if (ntr.eq.1) &
      nscm = 0

   nsolut =  nrp
   if ( nscm > 0 .and. ntb == 0 ) then
      ndfmin = 6   ! both translation and rotation com motion removed
      if (nsolut == 1) ndfmin = 3
      if (nsolut == 2) ndfmin = 5
   else if ( nscm > 0 ) then
      ndfmin = 3    ! just translation com will be removed
   else
      ndfmin = 0
   end if
   if (ibelly > 0) then   ! No COM Motion Removal, ever.
      nscm = 0
      ndfmin = 0

      !Do not allow ntt=3 with ibelly=1 - this does not make sense and will cause
      !issues.
      if (ntt==3) then
        call sander_bomb("mdread2","ibelly=1 with ntt=3 is not a valid option.", &
                         "Either use a different thermostat or avoid using ibelly.")
      end if
   end if
   if(nscm <= 0) nscm = 0
   if (therm_par > 0.0d0) ndfmin = 0 ! No COM motion removal for middle scheme NVT simulation
   if(gamma_ln > 0.0d0)ndfmin=0  ! No COM motion removal for LD simulation
   if(ntt == 4)ndfmin=0  ! No COM motion removal for Nose'-Hoover simulation
   init = 3
   if (irest > 0) init = 4
   if (dielc <= ZERO ) dielc = ONE
   if (tautp <= ZERO ) tautp = 0.2d0
   if (taup <= ZERO ) taup = 0.2d0

   !     ----- RESET THE CAP IF NEEDED -----

   ! ivcap == 0: Cap will be in effect if it is in the prmtop file (ifcap = 1)

   if(ivcap == 1) then
      ! Cap will be in effect even if not in prmtop file
      !   requires additional information in sander.in file as in the case of ivcap == 3, 4, or 5
      ifcap = 2
   else if(ivcap == 2) then
      ! Inactivate cap
      ifcap = 0
   else if(ivcap == 3) then
      ! Sphere -> not yet implemented
      ifcap = 3
   else if(ivcap == 4) then
      ! Orthorhombus
      ifcap = 4
   else if(ivcap == 5) then
      ! Shell of waters around solute
      ifcap = 5
   end if

   !Support for random seed using time of day in microsecs
   if( ig==-1 ) then
     !Turn on NO_NTT3_SYNC when ig=-1. This means the code no
     !longer synchronized the random numbers between streams when
     !running in parallel giving better scaling.
     no_ntt3_sync = 1
     ! call microsec(ig)
     call wallclock( wallc ) ! GNU fortran yields wallc to milliseconds
     wallc = modulo( 1.d3*wallc, 1.d6) ! should give six digits, positive
     ig = wallc
#ifdef MPI
     write (6, '(a,i8,a)') "Note: ig = -1. Setting random seed to ", ig ," based on wallclock &
                               &time in microseconds"
     write (6, '(a)') "      and disabling the synchronization of random &
                               &numbers between tasks"
     write (6, '(a)') "      to improve performance."
#else
#  ifndef API
     write (6, '(a,i8,a)') "Note: ig = -1. Setting random seed to ", ig ," based on wallclock &
                                &time in microseconds."
#  endif /* API */
#endif
   end if

#ifdef MPI
   !For some runs using multisander where the coordinates of the two 'replicas' need to be identical,
   !for example TI, it is critical that the random number stream is synchronized between all replicas.
   !Only use the IG value from worldrank=0. Ok to broadcast between the various sander masters since
   !they all call mdread2.
   !Also needs to be synchronized for adaptive QM/MM (qmmm_nml%vsolv > 1)
   if ( icfe > 0)  then
      ! no_ntt3_sync currently does not work with softcore TI simulations
      ! see sc_lngdyn in softcore.F90
      ! AWG: I think I also need to set no_ntt3_sync=0 ???
      if (ifsc > 0) no_ntt3_sync = 0
      call mpi_bcast(ig, 1, MPI_INTEGER, 0, commmaster, ierr)
   end if
#endif

   ! -------------------------------------------------------------------
   !     ----- PRINT DATA CHARACTERIZING THE RUN -----
   ! -------------------------------------------------------------------

   nr = nrp
#ifndef API
   write(6,9328)
   write(6,9008) title
   write(6,'(/a)') 'General flags:'
   write(6,'(5x,2(a,i8))') 'imin    =',imin,', nmropt  =',nmropt
#endif

   ! Error Checking for REMD
#ifdef MPI
   if (rem/=0) then
      ! Make sure that the number of replicas is even so
      !  that they all have partners in the exchange step
      if (mod(numgroup,2).ne.0) then
         write (6,'(a)') "==================================="
         write (6,'(a)') "REMD requires an even # of replicas"
         write (6,'(a)') "==================================="
         FATAL_ERROR
      endif

      write(6,'(/a)') 'Replica exchange'
      write(6,'(5x,4(a,i8))') 'numexchg=',numexchg,', rem=',rem

      ! REPCRD option temporarily disabled
      if (repcrd == 0) write(6,'(a)') &
         "REMD WARNING: repcrd disabled. Only replica &
         &trajectories/output can be written."

      ! Check for correct number of exchanges
      if (numexchg <= 0) then
         write(6,'(a)') "REMD ERROR: numexchg must be > 0, "
         FATAL_ERROR
      endif

      ! RXSGLD
      if (isgld > 0 .and. rem < 0) then
         write(6, '(a)') 'Multi-D REMD and replica-exchange SGLD are not &
                         &supported yet!'
         FATAL_ERROR
      end if

      ! Hybrid GB
      if (numwatkeep >= 0) then
         write(6,'(5x,4(a,i8))') 'numwatkeep=',numwatkeep,', hybridgb=',hybridgb
         ! Check that user specified GB model for hybrid REMD
         if (hybridgb /= 2 .and. hybridgb /= 5 .and. hybridgb /= 1) then
            write(6,'(a)') "HYBRID REMD ERROR: hybridgb must be 1, 2, or 5."
            FATAL_ERROR
         endif
      else
         !Check that user did not specify GB model if no hybrid run.
         if (hybridgb /= 0) then
            write(6,'(a)') &
            "HYBRID REMD ERROR: numwatkeep must be >= 0 when hybridgb is set."
            FATAL_ERROR
         endif
      endif
      ! RREMD
      if (rremd>0) then
         write(6,'(5x,4(a,i8))') "rremd=",rremd, ", reservoir_exchange_step=", &
                 & reservoir_exchange_step
         if (mod(reservoir_exchange_step, 2) .ne. 0) then
            write(6,'(a)') "reservoir_exchange_step must be a multiple of 2"
            FATAL_ERROR
         end if
      endif

      ! M-REMD (rem < 0) requires netcdf output.
      if (rem < 0 .and. ioutfm .ne. 1) then
         write(6,'(a)') "ERROR: Multi-D REMD (rem < 0) requires NetCDF &
                        &trajectories (ioutfm=1)"
         FATAL_ERROR
      endif

#  ifdef LES
      ! DAN ROE: Temporarily disable LES REMD until it is verified with new
      !           REMD code
      if (rem==2) then
         write (6,*) "******* LES REM (rem==2) temporarily disabled. Stop. *******"
         FATAL_ERROR
      endif

      if (rem==2 .and. igb/=1) then
         write (6,*) ' partial REM (rem=2) only works with igb=1'
         FATAL_ERROR
      endif
#  else
      if (rem==2) then
         write(6,*) '*******  For rem ==2, partial REM'
         write(6,*) 'use sander.LES with topology created by addles'
         FATAL_ERROR
      endif
#  endif

   endif ! rem>0

#else /* _NO_ MPI */
   ! Check if user set numexchg with no MPI
   if (numexchg>0) write(6,'(a)') &
      "WARNING: numexchg > 0 - for REMD run please recompile sander for &
      &parallel runs."

   ! Check if user set numwatkeep or hybridgb with no MPI - not sensible.
   if (numwatkeep>=0) write(6,'(a)') &
      "WARNING: numwatkeep >= 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."

   if (hybridgb>0) write(6,'(a)') &
      "WARNING: hybridgb > 0 - for hybrid REMD run please recompile &
      &sander for parallel runs."
#endif /* MPI */
   ! End error checking for REMD

#ifndef API
   write(6,'(/a)') 'Nature and format of input:'
   write(6,'(5x,4(a,i8))') 'ntx     =',ntx,', irest   =',irest, &
         ', ntrx    =',ntrx

   write(6,'(/a)') 'Nature and format of output:'
   write(6,'(5x,4(a,i8))') 'ntxo    =',ntxo,', ntpr    =',ntpr, &
         ', ntrx    =',ntrx,', ntwr    =',ntwr
   write(6,'(5x,3(a,i8))') ntwx    =',ntwx, &
         ', ntwv    =',ntwv,', ntwe    =',ntwe
   write(6,'(5x,2(a,i8),a,i7)') 'ioutfm  =',ioutfm, &
         ', ntwprt  =',ntwprt, &
         ', rbornstat=',rbornstat
   if (ntwf > 0) &
      write(6,'(5x, a,i8)') 'ntwf    =',ntwf
   write(6,'(/a)') 'Potential function:'
   write(6,'(5x,5(a,i8))') 'ntf     =',ntf,', ntb     =',ntb, &
         ', igb     =',igb,', nsnb    =',nsnb
   write(6,'(5x,3(a,i8))') 'ipol    =',ipol,', gbsa    =',gbsa, &
         ', iesp    =',iesp
   write(6,'(5x,3(a,f10.5))') 'dielc   =',dielc, &
         ', cut     =',cut,', intdiel =',intdiel
   if (lj1264 /= 0) &
      write(6, '(5x,a,i8)') 'lj1264  =',lj1264

   ! charge relocation
   if ( ifcr /= 0 ) then
      write(6,'(/a)') 'Charge relocation:'
      write(6,'(5x,2(a,i8))') 'cropt   =', cropt, &
                                ', crprintcharges=', crprintcharges
      write(6,'(5x,2(a,f10.5))') 'crcut   =', crcut, ', crskin  =', crskin
   end if

   if (( igb /= 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. igb /= 8) &
                                  ) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', gbalpha= ',gbalpha
      write(6,'(5x,3(a,f10.5))') 'gbbeta  =',gbbeta, &
            ', gbgamma =',gbgamma,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
   end if

   ! print output for igb=8
   if ( igb == 8 ) then
      write(6,'(5x,3(a,f10.5))') 'saltcon =',saltcon, &
            ', offset  =',offset,', surften =',surften
      write(6,'(5x,3(a,f10.5))') 'rdt     =',rdt, ', rgbmax  =',rgbmax, &
            '  extdiel =',extdiel
      write(6,'(5x,3(a,i8))') 'alpb  = ',alpb
      write(6,'(5x,3(a,f10.5))') 'gbalphaH  =',gbalphaH, &
            ', gbbetaH   =',gbbetaH,',  gbgammaH  = ',gbgammaH
      write(6,'(5x,3(a,f10.5))') 'gbalphaC  =',gbalphaC, &
            ', gbbetaC   =',gbbetaC,',  gbgammaC  = ',gbgammaC
      write(6,'(5x,3(a,f10.5))') 'gbalphaN  =',gbalphaN, &
            ', gbbetaN   =',gbbetaN,',  gbgammaN  = ',gbgammaN
      write(6,'(5x,3(a,f10.5))') 'gbalphaOS =',gbalphaOS, &
            ', gbbetaOS  =',gbbetaOS,',  gbgammaOS = ',gbgammaOS
      write(6,'(5x,3(a,f10.5))') 'gbalphaP  =',gbalphaP, &
            ', gbbetaP   =',gbbetaP,',  gbgammaP  = ',gbgammaP
      ! gbneck2nu pars
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_hnu  =', gb_alpha_hnu, &
        ', gb_beta_hnu   =', gb_beta_hnu, ',  gb_gamma_hnu  = ', gb_gamma_hnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_cnu  =', gb_alpha_cnu, &
        ', gb_beta_cnu   =', gb_beta_cnu, ',  gb_gamma_cnu  = ', gb_gamma_cnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_nnu  =', gb_alpha_nnu, &
        ', gb_beta_nnu   =', gb_beta_nnu, ',  gb_gamma_nnu  = ', gb_gamma_nnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_onu  =', gb_alpha_osnu, &
        ', gb_beta_onu   =', gb_beta_osnu, ',  gb_gamma_onu  = ', gb_gamma_osnu
      write(6,'(5x,3(a,f10.5))') 'gb_alpha_pnu  =', gb_alpha_pnu, &
        ', gb_beta_pnu   =', gb_beta_pnu, ',  gb_gamma_pnu  = ', gb_gamma_pnu
   end if


   if( alpb /= 0 ) then
      write(6,'(5x,3(a,f10.5))') 'Arad =', Arad
   end if

   write(6,'(/a)') 'Frozen or restrained atoms:'
   write(6,'(5x,4(a,i8))') 'ibelly  =',ibelly,', ntr     =',ntr
   if( ntr == 1 ) write(6,'(5x,a,f10.5)') 'restraint_wt =', restraint_wt

   if( imin /= 0 ) then

      write(6,'(/a)') 'Energy minimization:'
      ! print inputable variables applicable to all minimization methods.
      write(6,'(5x,4(a,i8))') 'maxcyc  =',maxcyc,', ncyc    =',ncyc, &
            ', ntmin   =',ntmin
      write(6,'(5x,2(a,f10.5))') 'dx0     =',dx0, ', drms    =',drms

      ! Input flag ntmin determines the method of minimization
      select case ( ntmin )
      case ( LMOD_NTMIN_XMIN, LMOD_NTMIN_LMOD )
         call write_lmod_namelist( )
      case default
         ! invalid ntmin
         write(6,'(/2x,a,i3,a)') 'Error: Invalid NTMIN (',ntmin,').'
         stop
      end select
   else
      write(6,'(/a)') 'Molecular dynamics:'
      write(6,'(5x,4(a,i10))') 'nstlim  =',nstlim,', nscm    =',nscm, &
            ', nrespa  =',nrespa
      write(6,'(5x,3(a,f10.5))') 't       =',t, &
            ', dt      =',dt,', vlimit  =',vlimit

      if( ntp /= 0 ) then
         write(6,'(/a)') 'Pressure regulation:'
         write(6,'(5x,4(a,i8))') 'ntp     =',ntp
         write(6,'(5x,3(a,f10.5))') 'pres0   =',pres0, &
               ', comp    =',comp,', taup    =',taup
         if (barostat == 2) then
            write(6, '(5x,a)') 'Monte-Carlo Barostat:'
            write(6, '(5x,a,i8)') 'mcbarint  =', mcbarint
         end if
      end if

      if (baroscalingdir /= 0) then
         write(6,'(/a)') 'Directional pressure scaling:'
         write(6,'(5x,a,i8)') 'baroscalingdir  =', baroscalingdir
      end if

   end if

   if( ntc /= 1 ) then
      write(6,'(/a)') 'SHAKE:'
      write(6,'(5x,4(a,i8))') 'ntc     =',ntc,', jfastw  =',jfastw
      write(6,'(5x,3(a,f10.5))') 'tol     =',tol
   end if

   if( ifcap == 1 .or. ifcap == 2 .or. ifcap == 3 ) then
      write(6,'(/a)') 'Water cap:'
      write(6,'(5x,2(a,i8))') 'ivcap   =',ivcap,', natcap  =',natcap
      write(6,'(5x,2(a,f10.5))') 'fcap    =',fcap, ', cutcap  =',cutcap
      write(6,'(5x,3(a,f10.5))') 'xcap    =',xcap, ', ycap    =',ycap,    &
                                 ', zcap    =',zcap
   else if( ifcap == 4 ) then
      write(6,'(/a)') 'Orthorhombus:'
      write(6,'(5x,1(a,i8))')    'ivcap   =',ivcap
      write(6,'(5x,1(a,f10.5))') 'forth   =',forth
      write(6,'(5x,3(a,f10.5))') 'xlorth  =',xlorth,', ylorth  =',ylorth, &
                                 ', zlorth  =',zlorth
      write(6,'(5x,3(a,f10.5))') 'xorth   =',xorth, ', yorth   =',yorth,  &
                                 ', zorth   =',zorth
   else if( ifcap == 5 ) then
      write(6,'(/a)') 'Water shell:'
      write(6,'(5x,(a,i8,a,f10.5))') 'ivcap   =',ivcap,', cutcap  =',cutcap
   endif

   if( nmropt > 0 ) then
      write(6,'(/a)') 'NMR refinement options:'
      write(6,'(5x,4(a,i8))')'iscale  =',iscale,', noeskp  =',noeskp, &
            ', ipnlty  =',ipnlty,', mxsub   =',mxsub
      write(6,'(5x,3(a,f10.5))') 'scalm   =',scalm, &
            ', pencut  =',pencut,', tausw   =',tausw
   end if

   if( numextra > 0 ) then
      write(6,'(/a)') 'Extra-points options:'
      write(6,'(5x,4(a,i8))') 'frameon =',frameon, &
            ', chngmask=',chngmask
   end if

   if( ipol > 0 ) then
      write(6,'(/a)') 'Polarizable options:'
      write(6,'(5x,4(a,i8))') 'indmeth =',indmeth, &
            ', maxiter =',maxiter,', irstdip =',irstdip, &
            ', scaldip =',scaldip
      write(6,'(5x,3(a,f10.5))') &
            'diptau  =',diptau,', dipmass =',dipmass
     if ( ipol > 1 ) then
      write(6,'(5x,3(a,f10.5))') &
            'Default Thole coefficient = ',dipdamp
     end if
   end if

#  ifdef MPI /* SOFT CORE */
   if( icfe /= 0 .or. ifsc/=0) then
      write(6,'(/a)') 'Free energy options:'
      write(6,'(5x,3(a,i8))')       'icfe    =', icfe   , ', ifsc    =', ifsc,    ', klambda =', klambda
      write(6,'(5x,3(a,f8.4))')     'clambda =', clambda, ', scalpha =', scalpha, ', scbeta  =', scbeta
      write(6,'(5x,2(a,i8))')       'sceeorder =', sceeorder, ' dvdl_norest =', dvdl_norest
      write(6,'(5x,a,f8.4,a,i8)')   'dynlmb =', dynlmb, ' logdvdl =', logdvdl
   end if
   if ( ifmbar /= 0 ) then
      write (6,'(/a)') 'FEP MBAR options:'
      write(6,'(5x,a,i8,a,i8)') 'ifmbar  =', ifmbar, ',  bar_intervall = ', bar_intervall
      write(6,'(5x,3(a,f6.4))') 'bar_l_min =', bar_l_min, ',  bar_l_max =', bar_l_max, &
                                ',  bar_l_incr =', bar_l_incr
   end if
#  endif /* MPI */

   if (ilrt /= 0) then
      write (6,*)
      write (6,'(a,i4,a,i4)') ' Linear Response Theory: ilrt =', ilrt, ' lrt_interval =', lrt_interval
      write (6,*)
   end if

   if( itgtmd == 1 ) then
      write(6,'(/a)') 'Targeted molecular dynamics:'
      write(6,'(5x,3(a,f10.5))') 'tgtrmsd =',tgtrmsd, &
            ', tgtmdfrc=',tgtmdfrc
   end if

   ntrelax = 0 ! needed for proper behavior of timing
   ntrelaxe = 0 ! needed for proper behavior of timing

   if( ntb > 0 ) then
      write(6,'(/a)') 'Ewald parameters:'
      write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
            ', ew_type =',ew_type,', nbflag  =',nbflag, &
            ', use_pme =',use_pme
      write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
            ', eedmeth =',eedmeth,', netfrc  =',netfrc
      write(6, 9002) a, b, c
      write(6, 9003) alpha, beta, gamma
      write(6, 9004) nfft1, nfft2, nfft3
      write(6, 9006) cutoffnb, dsum_tol
      write(6, 9007) ew_coeff
      write(6, 9005) order
      9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
      9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
      9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
      9005 format (5x,'Interpolation order =',i5)
      9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e10.3)
      9007 format (5x,'Ewald Coefficient =',f9.5)
   end if

#ifdef MPI
  ! Prints an error if REMD is to be performed with different cell sizes across
  ! the different replicas
  if (rem .ne. 0 .and. ntp .eq. 0 .and. ntb .ne. 0) then
    allocate(repvals(numgroup), stat=ierr)
    REQUIRE(ierr == 0)
    do i = 1, 6
       repvals(:) = 0.d0
       if (i .eq. 1) then
         val = a
       else if (i .eq. 2) then
         val = b
       else if (i .eq. 3) then
         val = c
       else if (i .eq. 4) then
         val = alpha
       else if (i .eq. 5) then
         val = beta
       else if (i .eq. 6) then
         val = gamma
       end if
       call mpi_allgather(val, 1, mpi_double_precision, &
                          repvals, 1, mpi_double_precision, &
                          commmaster, ierr)
       do j = 2, numgroup
         if (abs(repvals(j)-repvals(1)) .gt. 1.0d-5) then
           write(6,'(a,a,a)') error_hdr, "The box sizes do not match for all REMD replicas."
           write(6,'(11x,a)') "This can happen, for example, if you are using different input coordinate files"
           write(6,'(11x,a)') "for the different replicas and these structures come from a NPT run."
           write(6,'(11x,a)') "Make sure the box lengths and angles are the same for all replicas in all input"
           write(6,'(11x,a)') "coordinate files."
           FATAL_ERROR
         end if
       end do
    end do
  end if
#endif

!---- QMMM Options ----

   if( qmmm_nml%ifqnt ) then
      write(6, '(/a)') 'QMMM options:'
      write(6, '(5x,"        ifqnt = True       nquant = ",i8)') &
                 qmmm_struct%nquant
      write(6, '(5x,"         qmgb = ",i8,"  qmcharge = ",i8,"   adjust_q = ",i8)') &
                 qmmm_nml%qmgb, qmmm_nml%qmcharge, qmmm_nml%adjust_q
      write(6, '(5x,"         spin = ",i8,"     qmcut = ",f8.4, "    qmshake = ",i8)') qmmm_nml%spin, &
                 qmmm_nml%qmcut, qmmm_nml%qmshake
      write(6, '(5x,"     qmmm_int = ",i8)') qmmm_nml%qmmm_int
      write(6, '(5x,"lnk_atomic_no = ",i8,"   lnk_dis = ",f8.4," lnk_method = ",i8)') &
                 qmmm_nml%lnk_atomic_no,qmmm_nml%lnk_dis, qmmm_nml%lnk_method
      if ( qmmm_nml%qmtheory%PM3 ) then
         write(6, '(5x,"     qm_theory =     PM3")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%AM1 ) then
         write(6, '(5x,"     qm_theory =     AM1")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%AM1D ) then
         write(6, '(5x,"     qm_theory =   AM1/d")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%MNDO ) then
         write(6, '(5x,"     qm_theory =    MNDO")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%MNDOD ) then
         write(6, '(5x,"     qm_theory =  MNDO/d")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PDDGPM3 ) then
         write(6, '(5x,"     qm_theory = PDDGPM3")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PDDGMNDO ) then
         write(6, '(5x,"     qm_theory =PDDGMNDO")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PM3CARB1 ) then
         write(6, '(5x,"     qm_theory =PM3CARB1")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%DFTB ) then
         write(6, '(5x,"     qm_theory =    DFTB")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%RM1 ) then
         write(6, '(5x,"     qm_theory =     RM1")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PDDGPM3_08 ) then
         write(6, '(5x,"     qm_theory = PDDGPM3_08")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PM6 ) then
         write(6, '(5x,"     qm_theory =     PM6")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PM3ZNB ) then
         write(6, '(5x,"     qm_theory = PM3/ZnB")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%PM3MAIS ) then
         write(6, '(5x,"     qm_theory = PM3-MAIS")',ADVANCE='NO')
      else if ( qmmm_nml%qmtheory%EXTERN ) then
         write(6, '(5x,"     qm_theory =     EXTERN")',ADVANCE='NO')
      else
         write(6, '(5x,"     qm_theory = UNKNOWN!")',ADVANCE='NO')
      end if
      write (6, '(" verbosity = ",i8)') qmmm_nml%verbosity
      if (qmmm_nml%qmqm_analyt) then
        write(6, '(5x,"       qmqmdx = Analytical")')
      else
        write(6, '(5x,"       qmqmdx = Numerical")')
      end if
      !AWG: if EXTERN in use, skip printing of options that do not apply
      EXTERN: if ( .not. qmmm_nml%qmtheory%EXTERN ) then
         if (qmmm_nml%tight_p_conv) then
            write(6, '(5x," tight_p_conv = True (converge density to SCFCRT)")')
         else
            write(6, '(5x," tight_p_conv = False (converge density to 0.05xSqrt[SCFCRT])")')
         end if
         write(6, '(5x,"      scfconv = ",e10.3,"  itrmax = ",i8)') qmmm_nml%scfconv, qmmm_nml%itrmax
         if (qmmm_nml%printcharges) then
            write(6, '(5x," printcharges = True ")',ADVANCE='NO')
         else
            write(6, '(5x," printcharges = False")',ADVANCE='NO')
         end if
         select case (qmmm_nml%printdipole)
            case (1)
               write(6, '(5x," printdipole = QM   ")',ADVANCE='NO')
            case (2)
               write(6, '(5x," printdipole = QM+MM")',ADVANCE='NO')
            case default
               write(6, '(5x," printdipole = False")',ADVANCE='NO')
         end select
         if (qmmm_nml%peptide_corr) then
            write(6, '(5x," peptide_corr = True")')
         else
            write(6, '(5x," peptide_corr = False")')
         end if
         if (qmmm_nml%qmmmrij_incore) then
            write(6, '(4x,"qmmmrij_incore = True ")')
         else
            write(6, '(4x,"qmmmrij_incore = False")')
         end if
         if (qmmm_nml%qmqm_erep_incore) then
            write(6, '(2x,"qmqm_erep_incore = True ")')
         else
            write(6, '(2x,"qmqm_erep_incore = False")')
         end if
         if (qmmm_nml%allow_pseudo_diag) then
            write(6, '(7x,"pseudo_diag = True ")',ADVANCE='NO')
            write(6, '("pseudo_diag_criteria = ",f8.4)') qmmm_nml%pseudo_diag_criteria
         else
            write(6, '(7x,"pseudo_diag = False")')
         end if
         write(6,   '(6x,"diag_routine = ",i8)') qmmm_nml%diag_routine
      end if EXTERN
      !If ntb=0 or use_pme =0 then we can't do qm_ewald so overide what the user may
      !have put in the namelist and set the value to false.
      if (qmmm_nml%qm_ewald>0) then
        if (qmmm_nml%qm_pme) then
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = True ")') qmmm_nml%qm_ewald
        else
          write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
        end if
          write(6, '(10x,"  kmaxqx = ",i4," kmaxqy = ",i4," kmaxqz = ",i4," ksqmaxq = ",i4)') &
                  qmmm_nml%kmaxqx, qmmm_nml%kmaxqy, qmmm_nml%kmaxqz, qmmm_nml%ksqmaxq
      else
        write(6, '(10x,"qm_ewald = ",i8, " qm_pme = False ")') qmmm_nml%qm_ewald
      end if
      !Print the fock matrix prediction params if it is in use.
      if (qmmm_nml%fock_predict>0) then
          write(6, '(6x,"fock_predict = ",i4)') qmmm_nml%fock_predict
          write(6, '(6x,"    fockp_d1 = ",f8.4," fockp_d2 = ",f8.4)') qmmm_nml%fockp_d1, qmmm_nml%fockp_d2
          write(6, '(6x,"    fockp_d2 = ",f8.4," fockp_d4 = ",f8.4)') qmmm_nml%fockp_d3, qmmm_nml%fockp_d4
      end if

      if (qmmm_nml%qmmm_switch) then
         write(6, '(7x,"qmmm_switch = True",3x,"r_switch_lo =",f8.4,3x,"r_switch_hi =",f8.4)') &
               & qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi
      !else
      !   write(6, '(7x,"qmmm_switch = False")')
      end if

      if (qmmm_nml%printdipole==2) then
         write(6, '("|",2x,"INFO: To compute MM dipole WAT residues will be stripped")')
      end if
   end if

   if (qmmm_nml%vsolv > 0) then
      call print(qmmm_vsolv)
   end if

! ----EMAP Options-----
   if(temap)call emap_options(5)
! ---------------------

#  ifdef MPI
! --- MPI TIMING OPTIONS ---
      write(6, '(/a)') '| MPI Timing options:'
      write(6, '("|",5x," profile_mpi = ",i8)') profile_mpi
! Sanity check for profile_mpi
      call int_legal_range('profile_mpi',profile_mpi,0,1)
! --------------------------
#  endif /* MPI */
#endif /* API */

   cut = cut*cut
   cut_inner = cut_inner*cut_inner

   !------------------------------------------------------------------------
   ! If user has requested generalized born electrostatics, set up variables
   !------------------------------------------------------------------------

   if( igb == 0 .and. gbsa > 0 ) then
      write(0,*) 'GB/SA calculation is performed only when igb>0'
      FATAL_ERROR
   end if
   if( gbsa == 2 .and. &
       ((imin == 0 .and. nstlim > 1) .or. &
        (imin == 1 .and. maxcyc > 1)) ) then
      write(0,*) 'GBSA=2 only works for single point energy calc'
      FATAL_ERROR
   end if
   if( gbsa == 3 ) then
      write(0,*) 'GBSA=3 only works for pmemd, not sander'
      FATAL_ERROR
   end if
   if ( igb /= 0 .and. igb /= 10 .and. ipb == 0 ) then
#if defined (LES) && !defined (API)
      write(6,*) 'igb=1,5,7 are working with LES, no SA term included'
#endif
      ! igb7 uses special S_x screening params.
      ! overwrite the tinker values read from the prmtop
      if (igb == 7) then
         do i=1,natom
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            if (atomicnumber .eq. 6) then
               x(l96+i-1) = 4.84353823306d-1
            else if (atomicnumber .eq. 1) then
               x(l96+i-1) = 1.09085413633d0
            else if (atomicnumber .eq. 7) then
               x(l96+i-1) = 7.00147318409d-1
            else if (atomicnumber .eq. 8) then
               x(l96+i-1) = 1.06557401132d0
            else if (atomicnumber .eq. 16) then
               x(l96+i-1) = 6.02256336067d-1
            else if (atomicnumber .eq. 15) then
               x(l96+i-1) = 5d-1
            else
               x(l96+i-1) = 5d-1
            end if
         end do
      end if

      ! changing S_x screening params for igb = 8
      ! overwrite the tinker values read from the prmtop

      if (igb == 8) then
         do i=1,natom
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if

            call isnucat(nucat,i,nres,60,ix(i02:i02+nres-1),ih(m02:m02+nres-1))
            ! check if atom belong to nucleic or protein residue
            !60 = the total number of different nucleic acid residue type in AMBER
            !Need to update this number if adding new residue type
            !update isnucat function in gb_ene.F90 too
            if (nucat == 1) then
                ! If atom belongs to nuc, use gbneck2nu pars
                if (atomicnumber .eq. 6) then
                  x(l96+i-1) = screen_cnu
                else if (atomicnumber .eq. 1) then
                  x(l96+i-1) = screen_hnu
                else if (atomicnumber .eq. 7) then
                  x(l96+i-1) = screen_nnu
                else if (atomicnumber .eq. 8) then
                  x(l96+i-1) = screen_onu
                else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = screen_pnu
                else
                  x(l96+i-1) = 5.d-1 ! not optimized
                end if
            else
               ! protein
               if (atomicnumber .eq. 6) then
                  x(l96+i-1) = Sc
               else if (atomicnumber .eq. 1) then
                  x(l96+i-1) = Sh
               else if (atomicnumber .eq. 7) then
                  x(l96+i-1) = Sn
               else if (atomicnumber .eq. 8) then
                  x(l96+i-1) = So
               else if (atomicnumber .eq. 16) then
                  x(l96+i-1) = Ss
               else if (atomicnumber .eq. 15) then
                  x(l96+i-1) = Sp ! We still don't have an optimized Sp parameter
                                  ! if nuc, use gbneck2nu pars
               else
                  !for atom type Cl,Br,...
                  !These parameters are also not optimized.
                  x(l96+i-1) = 5d-1
               end if
            end if !check isnucat
         end do !loop i=1, natom
      end if ! igb = 8

      ! Set up for igb == 2, 5, 7, 8
      ! Put gb parameters in arrays
      if ( igb == 2 .or. igb == 5 .or. igb == 7 .or. &
           hybridgb == 2 .or. hybridgb == 5) then
        do i=1,natom
            x(l2402+i-1) = gbalpha
            x(l2403+i-1) = gbbeta
            x(l2404+i-1) = gbgamma
        end do
      end if

      ! IGB = 8, update gb_alpha, gb_beta, gb_gamma
       if ( igb == 8 ) then
           do i=1,natom
              if (ix(i100) .eq. 1) then
                 atomicnumber = ix(i100+i)
              else
                 call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                    atomicnumber, errFlag)
#ifdef LES
                 if(errFlag) then
                    call get_atomic_number(ih(m04+i-1), &
                       x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                       atomicnumber, errFlag)
                 end if
#endif
              end if

          ! check if atom belongs to protein or nucleic acid
          call isnucat(nucat,i,nres,60,ix(i02:i02+nres-1),ih(m02:m02+nres-1))
          if (nucat == 1) then
             !if atom belong to nucleic part, use nuc pars
             if (atomicnumber .eq. 1) then
               x(l2402+i-1) = gb_alpha_hnu
               x(l2403+i-1) = gb_beta_hnu
               x(l2404+i-1) = gb_gamma_hnu
             else if (atomicnumber .eq. 6) then
               x(l2402+i-1) = gb_alpha_cnu
               x(l2403+i-1) = gb_beta_cnu
               x(l2404+i-1) = gb_gamma_cnu
             else if (atomicnumber .eq. 7) then
               x(l2402+i-1) = gb_alpha_nnu
               x(l2403+i-1) = gb_beta_nnu
               x(l2404+i-1) = gb_gamma_nnu
             else if (atomicnumber .eq. 8) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 16) then
               x(l2402+i-1) = gb_alpha_osnu
               x(l2403+i-1) = gb_beta_osnu
               x(l2404+i-1) = gb_gamma_osnu
             else if (atomicnumber .eq. 15) then
               x(l2402+i-1) = gb_alpha_pnu
               x(l2403+i-1) = gb_beta_pnu
               x(l2404+i-1) = gb_gamma_pnu
             else
               ! Use GB^OBC (II) (igb = 5) set for other atoms
               x(l2402+i-1) = 1.0d0
               x(l2403+i-1) = 0.8d0
               x(l2404+i-1) = 4.851d0
             end if ! assign nucleic acid pars

          else
             !not nucleic part, use protein pars
             if (atomicnumber .eq. 1) then
                     x(l2402+i-1) = gbalphaH
                     x(l2403+i-1) = gbbetaH
                     x(l2404+i-1) = gbgammaH
             else if (atomicnumber .eq. 6) then
                     x(l2402+i-1) = gbalphaC
                     x(l2403+i-1) = gbbetaC
                     x(l2404+i-1) = gbgammaC
             else if (atomicnumber .eq. 7) then
                     x(l2402+i-1) = gbalphaN
                     x(l2403+i-1) = gbbetaN
                     x(l2404+i-1) = gbgammaN
             else if (atomicnumber .eq. 8) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 16) then
                     x(l2402+i-1) = gbalphaOS
                     x(l2403+i-1) = gbbetaOS
                     x(l2404+i-1) = gbgammaOS
             else if (atomicnumber .eq. 15) then
                     x(l2402+i-1) = gbalphaP
                     x(l2403+i-1) = gbbetaP
                     x(l2404+i-1) = gbgammaP
             else
                     !use GBOBC set for other atom types
                     x(l2402+i-1) = 1.0d0
                     x(l2403+i-1) = 0.8d0
                     x(l2404+i-1) = 4.851d0
             end if ! assign protein pars
          end if !check nucat
        end do !loop i=1,natom
       end if ! end IGB = 8 section

     !       put fs(i)*(rborn(i) - offset) into the "fs" array
      fsmax = 0.d0
      do i=1,natom
         x(l96-1+i) = x(l96-1+i)*( x(l97-1+i) - offset )
         fsmax = max( fsmax, x(l96-1+i) )
         if (rbornstat == 1) then
            x(l186-1+i) = 0.d0
            x(l187-1+i) = 999.d0
            x(l188-1+i) = 0.d0
            x(l189-1+i) = 0.d0
         end if
      end do

      !     ---------------------------------------------------------------------
      !       ---get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !         T = 298.15, epsext=78.5,

      kappa = sqrt( 0.10806d0 * saltcon )

      !       ---scale kappa by 0.73 to account(?) for lack of ion exclusions:

      kappa = 0.73d0* kappa

      !Set kappa for qmmm if needed
      qm_gb%kappa = kappa
      !     ---------------------------------------------------------------------

      if ( gbsa == 1 ) then

         !     --- assign parameters for calculating SASA according to the
         !         LCPO method ---

         do i=1,natom
            ix(i80+i-1)=0
         end do

         !         --- get the number of bonded neighbors for each atom:

         do i=1,nbona
            atom1=ix(iiba+i-1)/3+1
            atom2=ix(ijba+i-1)/3+1
            ix(i80+atom1-1)=ix(i80+atom1-1)+1
            ix(i80+atom2-1)=ix(i80+atom2-1)+1
         end do

         do i=1,natom
            ix(i80-i)=ix(i80+i-1)
         end do

         do i=1,nbonh
            atom1=ix(iibh+i-1)/3+1
            atom2=ix(ijbh+i-1)/3+1
            ix(i80-atom1)=ix(i80-atom1)+1
            ix(i80-atom2)=ix(i80-atom2)+1
         end do


         !         --- construct parameters for SA calculation; note that the
         !             radii stored in L165 are augmented by 1.4 Ang.

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            call upper(atype)
            if (ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            hybridization = ix(i80-i)
            nbond=ix(i80+i-1)
            if (atomicnumber .eq. 6) then
               if (hybridization .eq. 4) then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.56482d0
                     x(l175-1+i) = -0.19608d0
                     x(l180-1+i) = -0.0010219d0
                     x(l185-1+i) = 0.0002658d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.23348d0
                     x(l175-1+i) = -0.072627d0
                     x(l180-1+i) = -0.00020079d0
                     x(l185-1+i) = 0.00007967d0
                  else if (nbond == 4) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.00000d0
                     x(l175-1+i) = 0.00000d0
                     x(l180-1+i) = 0.00000d0
                     x(l185-1+i) = 0.00000d0
                  else
                     write(6,*) 'Unusual nbond for CT:', i, nbond, &
                        ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               else
                  if (nbond == 2) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.51245d0
                     x(l175-1+i) = -0.15966d0
                     x(l180-1+i) = -0.00019781d0
                     x(l185-1+i) = 0.00016392d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.070344d0
                     x(l175-1+i) = -0.019015d0
                     x(l180-1+i) = -0.000022009d0
                     x(l185-1+i) = 0.000016875d0
                  else
                     write(6,*) 'Unusual nbond for C :', i, nbond, &
                       ' Using default carbon LCPO parameters'
                     x(l165-1+i) = 1.70d0 + 1.4d0
                     x(l170-1+i) = 0.77887d0
                     x(l175-1+i) = -0.28063d0
                     x(l180-1+i) = -0.0012968d0
                     x(l185-1+i) = 0.00039328d0
                  end if
               end if
            else if (atomicnumber .eq. 8) then
               if (atype == 'O ') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.68563d0
                  x(l175-1+i) = -0.1868d0
                  x(l180-1+i) = -0.00135573d0
                  x(l185-1+i) = 0.00023743d0
               else if (atype == 'O2') then
                  x(l165-1+i) = 1.60d0 + 1.4d0
                  x(l170-1+i) = 0.88857d0
                  x(l175-1+i) = -0.33421d0
                  x(l180-1+i) = -0.0018683d0
                  x(l185-1+i) = 0.00049372d0
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.49392d0
                     x(l175-1+i) = -0.16038d0
                     x(l180-1+i) = -0.00015512d0
                     x(l185-1+i) = 0.00016453d0
                  else
                     write(6,*) 'Unusual nbond for O:', i, nbond, &
                        ' Using default oxygen LCPO parameters'
                     x(l165-1+i) = 1.60d0 + 1.4d0
                     x(l170-1+i) = 0.77914d0
                     x(l175-1+i) = -0.25262d0
                     x(l180-1+i) = -0.0016056d0
                     x(l185-1+i) = 0.00035071d0
                  end if
               end if
            else if(atomicnumber .eq. 7) then
               if (atype == 'N3') then
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.22599d0
                     x(l175-1+i) = -0.036648d0
                     x(l180-1+i) = -0.0012297d0
                     x(l185-1+i) = 0.000080038d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.051481d0
                     x(l175-1+i) = -0.012603d0
                     x(l180-1+i) = -0.00032006d0
                     x(l185-1+i) = 0.000024774d0
                  else
                     write(6,*) 'Unusual nbond for N3:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               else
                  if (nbond == 1) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.73511d0
                     x(l175-1+i) = -0.22116d0
                     x(l180-1+i) = -0.00089148d0
                     x(l185-1+i) = 0.0002523d0
                  else if (nbond == 2) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.41102d0
                     x(l175-1+i) = -0.12254d0
                     x(l180-1+i) = -0.000075448d0
                     x(l185-1+i) = 0.00011804d0
                  else if (nbond == 3) then
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.062577d0
                     x(l175-1+i) = -0.017874d0
                     x(l180-1+i) = -0.00008312d0
                     x(l185-1+i) = 0.000019849d0
                  else
                     write(6,*) 'Unusual nbond for N:', i, nbond, &
                        ' Using default nitrogen LCPO parameters'
                     x(l165-1+i) = 1.65d0 + 1.4d0
                     x(l170-1+i) = 0.078602d0
                     x(l175-1+i) = -0.29198d0
                     x(l180-1+i) = -0.0006537d0
                     x(l185-1+i) = 0.00036247d0
                  end if
               end if
            else if(atomicnumber .eq. 16) then
               if (atype == 'SH') then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.7722d0
                  x(l175-1+i) = -0.26393d0
                  x(l180-1+i) = 0.0010629d0
                  x(l185-1+i) = 0.0002179d0
               else
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.54581d0
                  x(l175-1+i) = -0.19477d0
                  x(l180-1+i) = -0.0012873d0
                  x(l185-1+i) = 0.00029247d0
               end if
            else if (atomicnumber .eq. 15) then
               if (nbond == 3) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               else if (nbond == 4) then
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.03873d0
                  x(l175-1+i) = -0.0089339d0
                  x(l180-1+i) = 0.0000083582d0
                  x(l185-1+i) = 0.0000030381d0
               else
                  write(6,*) 'Unusual nbond for P:', i, nbond, &
                     ' Using default phosphorus LCPO parameters'
                  x(l165-1+i) = 1.90d0 + 1.4d0
                  x(l170-1+i) = 0.3865d0
                  x(l175-1+i) = -0.18249d0
                  x(l180-1+i) = -0.0036598d0
                  x(l185-1+i) = 0.0004264d0
               end if
            else if (atype(1:1) == 'Z') then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atomicnumber .eq. 1) then
               x(l165-1+i) = 0.00000d0 + 1.4d0
               x(l170-1+i) = 0.00000d0
               x(l175-1+i) = 0.00000d0
               x(l180-1+i) = 0.00000d0
               x(l185-1+i) = 0.00000d0
            else if (atype == 'MG') then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(l165-1+i) = 1.18d0 + 1.4d0
               !  The following values were taken from O.sp3 with two bonded
               !  neighbors -> O has the smallest van der Waals radius
               ! compared to all other elements which had been parametrized
               x(l170-1+i) = 0.49392d0
               x(l175-1+i) = -0.16038d0
               x(l180-1+i) = -0.00015512d0
               x(l185-1+i) = 0.00016453d0
            else if (atype == 'F') then
               x(l165-1+i) = 1.47d0 + 1.4d0
               x(l170-1+i) = 0.68563d0
               x(l175-1+i) = -0.1868d0
               x(l180-1+i) = -0.00135573d0
               x(l185-1+i) = 0.00023743d0
            else
               ! write( 0,* ) 'bad atom type: ',atype
               ! call mexit( 6,1 )
               x(l165-1+i) = 1.70 + 1.4;
               x(l170-1+i) = 0.51245;
               x(l175-1+i) = -0.15966;
               x(l180-1+i) = -0.00019781;
               x(l185-1+i) = 0.00016392;
               write(6,'(a,a)') 'Using carbon SA parms for atom type', atype
            end if
         end do  !  i=1,natom
         !
      else if ( gbsa == 2 ) then

         !     --- assign parameters for calculating SASA according to the
         !         ICOSA method; the radii are augmented by 1.4 A ---

         do i=1,natom
            write(atype,'(a2)') ih(m06+i-1)
            if(ix(i100) .eq. 1) then
               atomicnumber = ix(i100+i)
            else
               call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
                  atomicnumber, errFlag)
#ifdef LES
               if(errFlag) then
                  call get_atomic_number(ih(m04+i-1), &
                     x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                     atomicnumber, errFlag)
               end if
#endif
            end if
            if (atomicnumber .eq. 7) then
               x(L165-1+i) = 1.55d0 + 1.4d0
            else if (atomicnumber .eq. 6) then
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 1) then
               x(L165-1+i) = 1.20d0 + 1.4d0
            else if (atomicnumber .eq. 8) then
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 15) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 16) then
               x(L165-1+i) = 1.80d0 + 1.4d0
            else if (atomicnumber .eq. 17) then
               ! Cl radius
               x(L165-1+i) = 1.70d0 + 1.4d0
            else if (atomicnumber .eq. 9) then
               ! F radius
               x(L165-1+i) = 1.50d0 + 1.4d0
            else if (atomicnumber .eq. 35) then
               ! Br radius
               ! Bondi, J. Phys. Chem. 1964, 68, 441.
               x(L165-1+i) = 1.85d0 + 1.4d0
            else if (atomicnumber .eq. 20) then
               ! Ca radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.33d0 + 1.4d0
            else if (atomicnumber .eq. 11) then
               ! Na radius
               ! Calculated from Aqvist, J. Phys. Chem. 1990, 94, 8021.
               x(L165-1+i) = 1.87d0 + 1.4d0
            else if (atomicnumber .eq. 30) then
               ! Zn radius
               ! Hoops, Anderson, Merz, JACS 1991, 113, 8262.
               x(L165-1+i) = 1.10d0 + 1.4d0
            else if (atomicnumber .eq. 12) then
               !  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
               !  Mg radius = 1.45A: Aqvist 1992
               x(L165-1+i) = 1.18d0 + 1.4d0
            else
               write( 0,* ) 'bad atom type: ',atype
               FATAL_ERROR
            end if

            ! dummy LCPO values:
            x(L170-1+i) = 0.0d0
            x(L175-1+i) = 0.0d0
            x(L180-1+i) = 0.0d0
            x(L185-1+i) = 0.0d0
            !  write(6,*) i,' ',atype,x(L165-1+i)
         end do  !  i=1,natom

      end if ! ( gbsa == 1 )

   end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )

   !-----------------------------------
   ! If a LRT calculation is requested,
   ! setup the icosa-SASA parameters
   ! this code is copied from above
   !-----------------------------------

   if ( ilrt /= 0 ) then
      !     --- assign parameters for calculating SASA according to the
      !         ICOSA method; the radii are augmented by 1.4 A ---
      do i=1,natom
         write(atype,'(a2)') ih(m06+i-1)
         if(ix(i100) .eq. 1) then
            atomicnumber = ix(i100+i)
         else
            call get_atomic_number(ih(m04+i-1), x(lmass+i-1), &
               atomicnumber, errFlag)
#ifdef LES
            if(errFlag) then
               call get_atomic_number(ih(m04+i-1), &
                  x(lmass+i-1)*lesfac((lestyp(i)-1)*nlesty+lestyp(i)), &
                  atomicnumber, errFlag)
            end if
#endif
         end if
         if (atomicnumber .eq. 7) then
            x(L165-1+i) = 1.55d0 + 1.4d0
         else if (atomicnumber .eq. 6) then
            x(L165-1+i) = 1.70d0 + 1.4d0
         else if (atomicnumber .eq. 1 .or. &
              ! added for lone pairs
              atype == 'EP') then
            x(L165-1+i) = 1.20d0 + 1.4d0
         else if (atomicnumber .eq. 8) then
            x(L165-1+i) = 1.50d0 + 1.4d0
         else if (atomicnumber .eq. 15) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 16) then
            x(L165-1+i) = 1.80d0 + 1.4d0
         else if (atomicnumber .eq. 12) then
            !             Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
            !             Mg radius = 1.45A: Aqvist 1992
            x(L165-1+i) = 1.18d0 + 1.4d0
         else
            write( 0,* ) 'bad atom type: ',atype,' cannot perform SASA calculation'
            FATAL_ERROR
         end if  !  atype(1:1) == 'N'
         x(L170-1+i) = 0.0d0
         x(L175-1+i) = 0.0d0
         x(L180-1+i) = 0.0d0
         x(L185-1+i) = 0.0d0
         !write(6,*) i,' ',atype,x(L165-1+i)
      end do  !  i=1,natom

   end if ! ( ilrt /= 0 )

   ! -------------------------------------------------------------------
   !
   ! -- add check to see if the space in nmr.h is likely to be
   !     too small for this run:
   !     [Note: this check does *not* indicate if MXTAU and MXP are
   !      too small.  There is no easy way to ensure this, since
   !      the experimental intensities are read in a namelist
   !      command: if too many intensities are input, the read
   !      statment may cause a coredump before returning control
   !      to the main program.  Be careful.  sigh....]

   if (natom > matom .and. nmropt > 1) then
      write(6,*) 'WARNING: MATOM in nmr.h is smaller than the ', &
            natom,' atoms in this molecule.'
      write(6,*) 'Printout of NMR violations may be compromised.'
   end if

   ! -------------------------------------------------------------------
   !     --- checks on bogus data ---
   ! -------------------------------------------------------------------

   inerr = 0

   if( icfe < 0 .or. icfe > 1 ) then
      write(6,*) 'icfe must be 0 or 1 (icfe=2 is no longer supported)'
      DELAYED_ERROR
   end if
   if( icfe /= 0 .and. numgroup /= 2 ) then
      write(6,*) 'numgroup must be 2 if icfe is set'
      DELAYED_ERROR
   end if
   if (ievb>0) then
#ifdef MPI
!KFW  if( numgroup /= 2 ) then
!KFW     write(6,*) 'numgroup must be 2 if ievb is set'
!KFW     DELAYED_ERROR
!KFW  end if
#else
      write(6,'(/2x,a)') 'Setting ievb>0 requires compilation with MPI'
      DELAYED_ERROR
#endif
   end if
   if( (igb > 0 .or. ipb /= 0) .and. numextra > 0) then
      if (igb /= 6) then
         write(6,'(a)') 'Cannot use igb>0 (except igb=6) with extra-point force fields'
         DELAYED_ERROR
      end if
   end if

!AMD validation
   if(iamd.gt.0)then
     if (EthreshD .eq. 0.d0 .and. alphaD .eq. 0.d0 .and. EthreshP .eq. 0.d0 .and. alphaP .eq. 0.d0) then
     write(6,'(a,i3)')'| AMD error all main parameters are 0.0 for Accelerated MD (AMD) or Windowed Accelerated MD (wAMD) '
     DELAYED_ERROR
     endif
    if(w_amd.gt.0)then
     if (EthreshD_w .eq. 0.d0 .and. alphaD_w .eq. 0.d0 .and. EthreshP_w .eq. 0.d0 .and. alphaP_w .eq. 0.d0) then
     write(6,'(a,i3)')'| AMD error all extra parameters are 0.0 for Windowed Accelerated MD (wAMD) LOWERING BARRIERS'
     DELAYED_ERROR
     endif
    endif
   endif

!GAMD not supported
   if (igamd .gt. 0) then
      write(6,*)'GaMD is not supported in Sander for now. Please use PMEMD instead.'
      inerr=1
   end if

   if (ips < 0 .or. ips > 6) then
      write(6,'(/2x,a,i3,a)') 'IPS (',ips,') must be between 0 and 6'
      DELAYED_ERROR
   end if
   if (ips /= 0 .and. ipol > 0 ) then
      write(6,'(/2x,a)') 'IPS and IPOL are inconsistent options'
      DELAYED_ERROR
   endif
   if (ips /= 0 .and. lj1264 > 0) then
      write(6, '(/2x,a)') 'IPS and the LJ 12-6-4 potential are incompatible'
      DELAYED_ERROR
   endif
   if ( (igb > 0 .or. ipb /= 0) .and. ips > 0 ) then
      write(6,'(/2x,a,i3,a,i3,a)') 'IGB (',igb,') and ips (',ips, &
          ') cannot both be turned on'
      DELAYED_ERROR
   end if
   if (igb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 &
         .and. igb /= 6 .and. igb /= 7 .and. igb /= 8 .and. igb /= 10) then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 0,1,2,5,6,7,8 or 10.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. alpb /= 1 )  then
      write(6,'(/2x,a,i3,a)') 'ALPB (',alpb,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (alpb /= 0 .and. igb /= 1 .and. igb /= 2 .and. igb /= 5 .and. igb /=7 )  then
      write(6,'(/2x,a,i3,a)') 'IGB (',igb,') must be 1,2,5, or 7 if ALPB > 0.'
      DELAYED_ERROR
   end if

   if (jar < 0 .or. jar > 1) then
      write(6,'(/2x,a,i3,a)') 'JAR (',jar,') must be 0 or 1'
      DELAYED_ERROR
   end if

#ifdef LES
   if( igb /= 0 .and. igb /= 1 .and. igb /= 5 .and. igb /=7 ) then
      write(6,'(/,a)') 'Error: LES is only compatible with IGB > 0,1,5,7'
      DELAYED_ERROR
   end if
   if( alpb /= 0) then
      write(6,'(/,a)') 'Error: LES is not compatible with ALPB'
      DELAYED_ERROR
   end if
   if( gbsa > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with GBSA > 0'
      DELAYED_ERROR
   end if
   if( qmmm_nml%ifqnt ) then
      write(6,'(/,a)') 'Error: LES is not compatible with QM/MM'
      DELAYED_ERROR
   end if
   if( ipol > 0 ) then
      write(6,'(/,a)') 'Error: LES is not compatible with IPOL > 0'
      DELAYED_ERROR
   end if
   if (temp0les >= 0.d0 .and. iscale > 0 ) then
      write (6,'(/,a)') 'Error: iscale cannot be used with temp0les'
      DELAYED_ERROR
   end if
#else
   rdt = 0
#endif
   if (irest /= 0 .and. irest /= 1) then
      write(6,'(/2x,a,i3,a)') 'IREST (',irest,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ibelly /= 0 .and. ibelly /= 1) then
      write(6,'(/2x,a,i3,a)') 'IBELLY (',ibelly,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (imin < 0) then
      write(6,'(/2x,a,i3,a)') 'IMIN (',imin,') must be >= 0.'
      DELAYED_ERROR
   end if
   if( imin>0 .and. ntmin>2 .and. iscale>0 ) then
      write(6, '(/,a)') 'cannot have iscale>0 when ntmin>2'
      DELAYED_ERROR
   end if
   if (ntxo .eq. 0) then
      write(6, '(/2x,a,a)') 'Old style binary restart files (ntxo=0) are no longer supported.'
      DELAYED_ERROR
   end if
   if (ntxo /= 1 .and. ntxo /= 2) then
      write(6, '(/2x,a,i3,a)') 'NTXO (',ntxo,') must be 1 or 2.'
      DELAYED_ERROR
   end if

   if (imin == 5) then
      if (ifbox /= 0 .and. ntb == 2) then
         write(6,'(/2x,a)') 'WARNING: IMIN=5 with changing periodic boundaries (NTB=2) can result in'
         write(6,'(/2x,a)') '         odd energies being calculated. Use with caution.'
      endif
#ifdef MPI
      if (sandersize > 1 .and. ntb == 2) then
         write(6,'(/2x,a)') 'ERROR: IMIN=5 and NTB=2 cannot be run with multiple processors.'
         DELAYED_ERROR
      endif
#endif
   end if


   if (fswitch > 0) then
      if (fswitch >= cut) then
         write(6, '(/2x,a)') 'fswitch must be less than the cutoff'
         DELAYED_ERROR
      end if
      if (igb /= 0 .or. ntb == 0) then
         write(6, '(/2x,a)') 'fswitch is only compatible with periodic simulations'
         DELAYED_ERROR
      end if
      if (lj1264 > 0) then
         write(6, '(/2x,a)') 'fswitch is incompatible with the 12-6-4 L-J model'
         DELAYED_ERROR
      end if
   end if

   if (iscale > mxvar) then
      write(6,9501) iscale,mxvar
      9501 format('ERROR: ISCALE (',i5,') exceeds MXVAR (',i5, &
            '). See nmr.h')
      DELAYED_ERROR
   end if
   if (ntx < 1 .or. ntx > 7) then
      write(6,'(/2x,a,i3,a)') 'NTX (',ntx,') must be in 1..7'
      DELAYED_ERROR
   end if

   if (ntb /= 0 .and. ntb /= 1 .and. ntb /= 2) then
      write(6,'(/2x,a,i3,a)') 'NTB (',ntb,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if
   if (ntt < 0 .or. ntt > 10) then                                      ! APJ
      write(6,'(/2x,a,i3,a)') 'NTT (',ntt,') must be between 0 and 10.' ! APJ
      DELAYED_ERROR
   end if
   if (ntt == 1 .and. tautp < dt) then
      write(6, '(/2x,a,f6.2,a)') 'TAUTP (',tautp,') < DT (step size)'
      DELAYED_ERROR
   end if
   if( ntt < 3 .or. ntt > 10) then  ! APJ
      if( gamma_ln > 0.d0 ) then
         write(6,'(a)') 'ntt must be 3 to 10 if gamma_ln > 0' ! APJ
         DELAYED_ERROR
      end if
   end if
   if (infe /= 0 .and. infe /= 1) then
      write(6,'(/2x,a,i3,a)') 'IFE (',infe,') must be 0 or 1.'
      DELAYED_ERROR
   end if

   if (ntp /= 0 .and. ntp /= 1 .and. ntp /= 2 .and. ntp /= 3) then
      write(6,'(/2x,a,i3,a)') 'NTP (',ntp,') must be 0, 1, 2, or 3.'
      DELAYED_ERROR
   end if
   if (ntp > 0 .and. taup < dt .and. barostat == 1) then
      write(6, '(/2x,a,f6.2,a)') 'TAUP (',taup,') < DT (step size)'
      DELAYED_ERROR
   end if
   if (npscal < 0 .or. npscal > 1) then
      write(6,'(/2x,a,i3,a)') 'NPSCAL (',npscal,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntp > 0) then
      if (barostat /= 1 .and. barostat /= 2) then
         write(6, '(/2x,a,i3,a)') 'BAROSTAT (', barostat, ') must be 1 or 2'
         DELAYED_ERROR
      end if
      if (barostat == 2) then
         if (mcbarint <= 0) then
            write(6, '(/2x,a,i3,a)') 'MCBARINT (',mcbarint,') must be positive'
            DELAYED_ERROR
         end if
         if (mcbarint >= nstlim) then
            write(6, '(a)') 'WARNING: mcbarint is greater than the number of &
                     &steps. This is effectively constant volume.'
         end if
      end if
   end if

   if (ntc < 1 .or. ntc > 4) then
      write(6,'(/2x,a,i3,a)') 'NTC (',ntc,') must be 1,2,3 or 4.'
      DELAYED_ERROR
   end if
   if (jfastw < 0 .or. jfastw > 4) then
      write(6,'(/2x,a,i3,a)') 'JFASTW (',jfastw,') must be 0->4.'
      DELAYED_ERROR
   end if

   if (ntf < 1 .or. ntf > 8) then
      write(6,'(/2x,a,i3,a)') 'NTF (',ntf,') must be in 1..8.'
      DELAYED_ERROR
   end if

   if (ioutfm /= 0 .and. ioutfm /= 1) then
      write(6,'(/2x,a,i3,a)') 'IOUTFM (',ioutfm,') must be 0 or 1.'
      DELAYED_ERROR
   end if

   if (ntpr < 0) then
      write(6,'(/2x,a,i3,a)') 'NTPR (',ntpr,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwx < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWX (',ntwx,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntwv < -1) then
      write(6,'(/2x,a,i3,a)') 'NTWV (',ntwv,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwf < -1) then
      write(6, '(/2x,a,i3,a)') 'NTWF (',ntwf,') must be >= -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ioutfm /= 1) then
      write (6, '(/2x,a)') 'IOUTFM must be 1 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwf == -1 .and. ioutfm /= 1) then
      write(6, '(/2x,a)') 'IOUTFM must be 1 for NTWF == -1.'
      DELAYED_ERROR
   end if
   if (ntwv == -1 .and. ntwx == 0) then
      write (6, '(/2x,a)') 'NTWX must be > 0 for NTWV == -1.'
      DELAYED_ERROR
   end if
   if (ntwe < 0) then
      write(6,'(/2x,a,i3,a)') 'NTWE (',ntwe,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntave < 0) then
      write(6,'(/2x,a,i3,a)') 'NTAVE (',ntave,') must be >= 0.'
      DELAYED_ERROR
   end if
   if (ntr /= 0 .and. ntr /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTR (',ntr,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (ntrx /= 0 .and. ntrx /= 1) then
      write(6,'(/2x,a,i3,a)') 'NTRX (',ntrx,') must be 0 or 1.'
      DELAYED_ERROR
   end if
   if (nmropt < 0 .or. nmropt > 2) then
      write(6,'(/2x,a,i3,a)') 'NMROPT (',nmropt,') must be in 0..2.'
      DELAYED_ERROR
   end if

   ! check settings related to ivcap

   if(ivcap == 3 .or. ivcap == 4) then
      write(6,'(/2x,a)') 'IVCAP == 3 and IVCAP == 4 currently not implemented'
      DELAYED_ERROR
   endif
   if (ivcap < 0 .and. ivcap > 5) then
      write(6,'(/2x,a)') 'IVCAP must be 0 ... 5'
      DELAYED_ERROR
   end if
   if ((ivcap == 1 .or. ivcap == 5) .and. igb /= 10 .and. ipb == 0) then
      write(6,'(/2x,a)') 'IVCAP == 1,5 only works with Poisson Boltzmann (igb=10 or ipb/=0)'
      DELAYED_ERROR
   end if
   if((ivcap == 1 .or. ivcap == 3 .or. ivcap == 5 ) .and. cutcap <= 0.0d0) then
      write(6,'(/2x,a)') 'For IVCAP == 1,3, or 5, cutcap must be > 0.0'
      DELAYED_ERROR
   endif
   if (ivcap == 4 .and. &
         (xlorth < ZERO .or. ylorth < ZERO .or. zlorth < ZERO .or. &
!  give magic numbers a name  srb aug 2007 !
          xorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          yorth > 47114710.0d0 .or. &
!  give magic numbers a name !
          zorth > 47114710.0d0)) then
      write(6,'(/2x,a)') &
      'For IVCAP == 4, xlorth, ylorth, zlorth, xorth, yorth, zorth must be set'
      DELAYED_ERROR
   end if
   if ((ivcap == 3 .or. ivcap == 4) .and. ibelly == 0) then
      write(6,'(/2x,a,a)') &
         'For IVCAP == 3 or 4, ibelly must be 1 and all atoms', &
         '  not in the spherical or orthorhombic region must be set NOT moving'
      DELAYED_ERROR
   end if
   if (ivcap == 5 .and. (imin /= 1 .or. maxcyc > 1)) then
      write(6,'(/2x,a,a)') &
         'IVCAP == 5 only works for single-point energy calculation'
      DELAYED_ERROR
   end if

   ! check if ifbox variable from prmtop file matches actual angles:

   if ( igb == 0 .and. ipb == 0 .and. ntb /= 0 ) then
      if ( ifbox == 1 ) then
         if ( abs(alpha - 90.0d0) > 1.d-5 .or. &
           abs(beta  - 90.0d0) > 1.d-5 .or. &
           abs(gamma - 90.0d0) > 1.d-5 ) then
           ifbox =3
#ifndef API
           write(6,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
#endif
         end if
      end if

      if ( ifbox == 2 ) then
         if ( abs(alpha - 109.4712190d0) > 1.d-5 .or. &
              abs(beta  - 109.4712190d0) > 1.d-5 .or. &
              abs(gamma - 109.4712190d0) > 1.d-5 ) then
              write(6,'(/2x,a)') &
              'Error: ifbox=2 in prmtop but angles are not correct'
              DELAYED_ERROR
         end if
      end if
   end if

   ! checks for targeted MD
   if (itgtmd /= 0 .and. itgtmd > 2) then
      write(6,'(/2x,a,i3,a)') 'ITGTMD (',itgtmd,') must be 0, 1 or 2.'
      DELAYED_ERROR
   end if
   if (itgtmd == 1 .and. ntr == 1) then
      if (len_trim(tgtfitmask) > 0 .or. len_trim(tgtrmsmask) <= 0) then
         write(6,'(/2x,a)') 'ITGTMD: tgtrmsmask (and not tgtfitmask) ' //  &
                            'should be specified if NTR=1'
         DELAYED_ERROR
      end if
   end if
   ! skip this test until fallback to rgroup() is supported
   !if (itgtmd == 1 .and. ntr == 0) then
   !   if (len_trim(tgtfitmask) == 0 .and. len_trim(tgtrmsmask) == 0) then
   !      write(6,'(/2x,a)')  &
   !        'ITGTMD: both tgtfitmask and tgtrmsmask should be specified if NTR=0'
   !      DELAYED_ERROR
   !   end if
   !end if

   !     -- consistency checking

   if (imin > 0.and.nrespa > 1)  then
      write(6,'(/2x,a)') 'For minimization, set nrespa,nrespai=1'
      DELAYED_ERROR
   end if
   if (ntp > 0 .and. nrespa > 1) then
      write(6,'(/2x,a)') 'nrespa must be 1 if ntp>0'
      DELAYED_ERROR
   end if
   if  (ntx < 4.and.init /= 3)  then
      write(6,'(/2x,a)') 'NTX / IREST inconsistency'
      DELAYED_ERROR
   end if
   if (ntb == 2 .and. ntp == 0) then
      write(6,'(/2x,a)') 'NTB set but no NTP option (must be 1 or 2)'
      DELAYED_ERROR
   end if
   if (ntp /= 0 .and. ntb /= 2) then
      write(6,'(/,a,a)')' NTP > 0 but not constant pressure P.B.C.', &
            ' (NTB = 2) must be used'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. ifbox == 0 .and. ntp /= 0) then
      write(6,'(/,a)') ' (NTB /= 0 and NTP /= 0) but IFBOX == 0'
      write(6,'(/,a)') ' This combination is not supported'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. &
         ( box(1) < 1.d0  .or. &
         box(2) < 1.d0  .or. &
         box(3) < 1.d0 ) ) then
      write(6,'(/,a,3f10.3)') ' BOX is too small: ',box(1),box(2),box(3)
      DELAYED_ERROR
   else if (ntb /= 0 .and. &
         (sqrt(cut) >= box(1)*0.5d0 .or. &
         sqrt(cut) >= box(2)*0.5d0 .or. &
         sqrt(cut) >= box(3)*0.5d0) ) then
      write(6,'(/,a)') ' CUT must be < half smallest box dimension'
      DELAYED_ERROR
   end if
   if (ntb /= 0 .and. (igb > 0 .or. ipb /= 0) ) then
      write(6,'(/,a)') ' igb>0 is only compatible with ntb=0'
      DELAYED_ERROR
   end if
   if ( ntb == 0 .and. sqrt(cut) < 8.05 .and. igb /= 10 .and. ipb == 0 .and. &
        igb /= 6 ) then
      write(6,'(/,a,f8.2)') ' unreasonably small cut for non-periodic run: ', &
         sqrt(cut)
      DELAYED_ERROR
   end if
   if ( rgbmax < 5.d0*fsmax ) then
      write(6,'(/,a,f8.2)') ' rgbmax must be at least ', 5.d0*fsmax
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. indmeth == 3 ) then
      write(6,'(/,a)') ' indmeth=3 cannot be used with icfe>0'
      DELAYED_ERROR
   end if
   if (icfe /= 0 .and. ibelly /= 0 ) then
      write(6,'(/,a)') ' ibelly cannot be used with icfe'
      DELAYED_ERROR
   end if
#ifdef MPI /* SOFT CORE */
   if (icfe /= 0 .and. dvdl_norest /= 0 ) then
      write(6,'(/,a)') 'dvdl_norest must == 0!'
      write(6,'(/,a)') 'The dvdl_norest option is deprecated.'
      write(6,'(/,a)') 'Restraint energies are seperated, &
        &and do not contribute to dvdl.'
      DELAYED_ERROR
   end if
#endif
   ! Modification done by Ilyas Yildirim
   if (icfe == 1 .and. (klambda < 1 .or. klambda > 6)) then
     write(6,'(/,a)') ' klambda must be between 1 and 6'
     DELAYED_ERROR
   end if
   ! End of modification done by Ilyas Yildirim

   if (clambda < 0.d0 .or. clambda > 1.d0 ) then
      write(6,'(/,a)') ' clambda must be between 0 and 1'
      DELAYED_ERROR
   end if

#ifdef MPI /* SOFT CORE */
   if (ifsc /= 0) then
      if (ifsc == 2) then
         write (6,'(/,a)') 'The ifsc=2 option is no longer supported. &
               &Internal energies of the soft core region.'
         write (6,'(a)') 'are now handled implicitly and setting ifsc=2 &
               &is no longer needed'
         DELAYED_ERROR
      end if
      if (icfe /= 1 .and. ifsc==1) then
         write (6,'(/,a)') ' Softcore potential requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if ( igb > 0 .or. ipb /= 0 ) then
         write (6,'(/,a)') ' Softcore potential is incompatible with GB (for now)'
         DELAYED_ERROR
      end if
      if ( ntf > 1 ) then
         write (6,'(/,a)') ' Softcore potentials require ntf=1 because SHAKE &
               &constraints on some bonds might be removed'
         DELAYED_ERROR
      end if
      if (clambda > 0.995 .or. clambda < 0.005) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with clambda < 0.005 or > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' Softcore potential requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
      if (imin == 1 .and. ntmin /= 2) then
         write (6,'(/,a)') ' Minimizations with ifsc=1 require the steepest descent algorithm.'
         write (6,'(/,a)') ' Set ntmin to 2 and restart'
         DELAYED_ERROR
      end if
   end if
   ! Following is temporary (dec 2014), until/unless problem gets fixed
   if ( logdvdl /= 0 .and. ifsc == 0 ) then
      write (6,'(/,a)') ' Cannot request logdvdl unless ifsc>0'
      DELAYED_ERROR
   end if
   if (ifmbar /= 0) then
      if (icfe /= 1) then
         write (6,'(/,a)') ' MBAR requires a standard TI run, set icfe to 1'
         DELAYED_ERROR
      end if
      if (ifsc /=0 .and. (bar_l_max > 0.995 .or. bar_l_min < 0.005) ) then
         write (6,'(/,a)') ' Softcore potentials cannot be used with &
               &bar_l_min < 0.005 or bar_l_max > 0.995'
         DELAYED_ERROR
      end if
      if (klambda /= 1) then
         write (6,'(/,a)') ' MBAR requires linear mixing, set klambda to 1'
         DELAYED_ERROR
      end if
   end if
#endif

   if (ilrt /= 0 .and. lrtmask == '') then
      write (6,'(a)') 'Linear Response Theory activated, but lrtmask is not set'
      DELAYED_ERROR
   end if

   !-----------------------------------------------------
   !     ----sanity checks for Ewald
   !-----------------------------------------------------

   if( igb == 0 .and. ipb == 0 ) then
      call float_legal_range('skinnb: (nonbond list skin) ', &
            skinnb,skinlo,skinhi)

      !  --- Will check on sanity of settings after coords are read in
      !      and the extent of the system is determined.

      if(periodic == 1)then
         call float_legal_range('skinnb+cutoffnb: (nonbond list cut) ', &
               skinnb+cutoffnb,zero,sphere)
      end if
      if (ntb==0 .and. use_pme/=0) then
         write(6,'(/,a)') &
         'Using PME with a non-periodic simulation does not make sense. Set either ntb>0 of use_pme=0.'
         DELAYED_ERROR
      end if
      call float_legal_range('a: (unit cell size) ',a,boxlo,boxhi)
      call float_legal_range('b: (unit cell size) ',b,boxlo,boxhi)
      call float_legal_range('c: (unit cell size) ',c,boxlo,boxhi)
      call float_legal_range('alpha: (unit cell angle) ', &
            alpha,anglo,anghi)
      call float_legal_range('beta: (unit cell angle) ', &
            beta,anglo,anghi)
      call float_legal_range('gamma: (unit cell angle) ', &
            gamma,anglo,anghi)
      call int_legal_range('order: (interpolation order) ', &
            order,orderlo,orderhi)
      call opt_legal_range('verbose: ',verbose,0,4)
      call opt_legal_range('netfrc: ',netfrc,0,1)
      call opt_legal_range('nbflag: ',nbflag,0,1)
      call opt_legal_range('nbtell: ',nbtell,0,2)
      call opt_legal_range('ew_type: ',ew_type,0,1)
      call opt_legal_range('vdwmeth: ',vdwmeth,0,2)
      call opt_legal_range('eedmeth: ',eedmeth,1,6)
      call opt_legal_range('ee_type: ',ee_type,1,2)
      call opt_legal_range('maxiter: ',maxiter,1,50)
      call opt_legal_range('indmeth: ',indmeth,0,3)
      call opt_legal_range('fix_quad: ',fix_quad,0,1)
      call float_legal_range('eedtbdns: (erfc table density) ', &
            eedtbdns,denslo,denshi)
   end if  ! ( igb == 0 .and. ipb == 0 )

   ! ---WARNINGS:

   if ( ibelly == 1 .and. igb == 0 .and. ipb == 0 .and. ntb /= 0 ) then
#ifndef API
      write(6,'(/,a,/,a,/,a)') &
            'Warning: Although EWALD will work with belly', &
            '(for equilibration), it is not strictly correct!'
#endif
   end if

   if (inerr == 1) then
      write(6, '(/,a)') ' *** input error(s)'
      FATAL_ERROR
   end if

   ! Load the restrained atoms (ntr=1) or the belly atoms (ibelly=1)
   ! or atoms for targeted md (itgtmd=1). Selections are read from
   ! &cntrl variables or, if these are not defined, it falls back to
   ! the old group input format.

   if(mtmd /= 'mtmd') then
     itgtmd = 2
     ntr = 0
     emtmd = 0.0d0
   end if
   konst = ntr > 0
   dotgtmd = itgtmd > 0
   belly = .false.
   natc = 0
   ngrp = 0
   natbel = 0
   nattgtfit = 0  ! number of atoms for tgtmd fitting (=overlap)
   nattgtrms = 0  ! number of atoms for tgtmd rmsd calculation
   nrc = 0
   if(konst.or.dotgtmd) then

      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph,solve)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#else
#  ifndef API
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.,solvph,solve)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#  endif /* API */
#endif /* LES */

      if(itgtmd == 2) then
        call mtmdcall(emtmd,x(lmtmd01),ix(imtmd02),x(lcrd),x(lforce),ih(m04),ih(m02),ix(i02),&
                    ih(m06),x(lmass),natom,nres,'READ')
      else
#ifndef API
        if (konst) write(6,9408)
        if (dotgtmd) write(6,9409)
#endif
        call rdrest(natom,ntrx,refc,x(lcrdr))
        !close(10)

        ! VH - tgtmd change: preferably call atommask() instead of rgroup()
        if (konst) then
#ifndef API
          if( len_trim(restraintmask) <= 0 ) then
              call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                          ih(m06),ih(m08),ix(icnstrgp),jgroup,indx,irespw, &
                          npdec,x(l60),konst,dotgtmd,belly,5,.true.)
          else
#endif
              if (mask_from_ref > 0) then ! Antoine Marion : base mask on reference coordinates
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                  ix(i02), ih(m02), x(lcrdr), restraintmask, ix(icnstrgp) )
              else
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                  ix(i02), ih(m02), x(lcrd), restraintmask, ix(icnstrgp) )
              endif

              ! for now, emulate the "GATHER ALL THE CONSTRAINED ATOMS TOGETHER"
              ! section of rgroup(); later, the various masks should be done
              ! differently, i.e. without the "gather", as in the following:
              !     x(l60:l60+natom-1) = restraint_wt
              !     natc = sum(ix(icnstrgp:icnstrgp+natom-1))

              natc = 0
              do i=1,natom
                if( ix(icnstrgp-1+i) <= 0 ) cycle
                natc = natc + 1
                ix(icnstrgp-1+natc) = i
                x(l60-1+natc) = restraint_wt
              end do
#ifndef API
              write(6,'(a,a,a,i5,a)') '     Mask ', &
              restraintmask(1:len_trim(restraintmask)), ' matches ',natc,' atoms'
          end if
#endif
        end if
        nrc = natc

        if (itgtmd == 1) then
#ifndef API
          if (len_trim(tgtfitmask) <= 0 .and. len_trim(tgtrmsmask) <= 0) then
              ! the following if-endif can be deleted when we stop
              ! supporting rgroup()
              if (konst) then
                ! cannot do both ntr and tgtmd together using old group format
                write(6,'(/2x,a)') 'NTR must be 0 for targeted MD (TGTMD=1)'
                call mexit(6,1)
              else  ! the following only for backward compatibility
                call rgroup(natom,natc,nres,ngrp,ix(i02),ih(m02),ih(m04), &
                            ih(m06),ih(m08),ix(icnstrgp),jgroup,indx,irespw, &
                            npdec,x(l60),konst,dotgtmd,belly,5,.true.)
                ! tgtmd atoms are now stored in nattgt, igroup -> icnstrgp
                nattgtfit = natc
                nattgtrms = natc
                do i=1,nattgtfit
                    ix(itgtfitgp-1+i) = ix(icnstrgp-1+i)
                    ix(itgtrmsgp-1+i) = ix(icnstrgp-1+i)
                end do
              end if
          else
#endif
              if (ntr == 0) then  ! read tgtfitmask only if ntr=1
                ! read in atom group for tgtmd fitting (=overlap region)
                call atommask( natom, nres, 0, ih(m04), ih(m06), &
                    ix(i02), ih(m02), x(lcrd), tgtfitmask, ix(itgtfitgp) )
                ! see comments above (for ntr) for the following reduction cycle
                nattgtfit = 0
                do i=1,natom
                  if( ix(itgtfitgp-1+i) <= 0 ) cycle
                  nattgtfit = nattgtfit + 1
                  ix(itgtfitgp-1+nattgtfit) = i
                end do
                write(6,'(a,a,a,i5,a)')  &
                '     Mask "', tgtfitmask(1:len_trim(tgtfitmask)-1),  &
                '" matches ',nattgtfit,' atoms'
              end if
              ! read in atom group for tgtmd rmsd calculation
              call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), tgtrmsmask, ix(itgtrmsgp) )
              nattgtrms = 0
              do i=1,natom
                if( ix(itgtrmsgp-1+i) <= 0 ) cycle
                nattgtrms = nattgtrms + 1
                ix(itgtrmsgp-1+nattgtrms) = i
              end do
              write(6,'(a,a,a,i5,a)')  &
              '     Mask "', tgtrmsmask(1:len_trim(tgtrmsmask)-1),  &
              '" matches ',nattgtrms,' atoms'
#ifndef API
          end if
#endif
        end if

      end if ! (itgtmd == 2)

   end if  ! (konst.or.dotgtmd)

   ! dotgtmd may be false here even if doing tgtmd
   ! this is so belly info is read properly? following existing KONST code

   dotgtmd=.false.
   konst = .false.
   belly = ibelly > 0
   ngrp = 0
   if(belly) then
      ! inserted here to fix the bug that coords are not available
      ! yet when distance based selection (<,>) is requested
#ifdef LES
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0les,.FALSE.,solvph,solve)
#  else
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#  endif
#else
#  ifdef MPI
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,temp0,.FALSE.,solvph,solve)
#  else
#    ifndef API
      call getcor(nr,x(lcrd),x(lvel),x(lforce),ntx,box,irest,t,.FALSE.)
#    endif
#  endif
#endif /* LES */
#ifndef API
      write(6,9418)
      if( len_trim(bellymask) <= 0 ) then
         call rgroup(natom,natbel,nres,ngrp,ix(i02),ih(m02),ih(m04),ih(m06), &
                     ih(m08),ix(ibellygp),jgroup,indx,irespw,npdec, &
                     x(l60),konst,dotgtmd,belly,5,.true.)
      else
#endif /* API */
         if (mask_from_ref > 0) then ! Antoine Marion : base mask on reference coordinates
           call rdrest(natom,ntrx,refc,x(lcrdr))
           call atommask( natom, nres, 0, ih(m04), ih(m06), &
              ix(i02), ih(m02), x(lcrdr), bellymask, ix(ibellygp) )
         else
           call atommask( natom, nres, 0, ih(m04), ih(m06), &
              ix(i02), ih(m02), x(lcrd), bellymask, ix(ibellygp) )
         endif
         natbel = sum(ix(ibellygp:ibellygp+natom-1))
#ifndef API
         write(6,'(a,a,a,i5,a)') '     Mask ', &
            bellymask(1:len_trim(bellymask)), ' matches ',natbel,' atoms'
      end if
#endif
   end if
   call setvar(ix,belly)

   !  see if the user has input a noshakemask string, and process it:
   natnos = 0
   if( len_trim(noshakemask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), noshakemask, noshakegp )
      natnos = sum(noshakegp(1:natom))
      write(6,*)
      write(6,'(a,a,a,i5,a)') 'Noshake mask ', &
         noshakemask(1:len_trim(noshakemask)), ' matches ',natnos,' atoms'
      call setnoshake(ix,noshakegp,ntc,num_noshake)
      if( ntf > 1 ) then
         write(6,'(a)') '   Setting ntf to 1'
         ntf = 1
      end if
   end if

#ifdef MPI /* SOFT CORE */
   ! lower charges if a crgmask is set
   if ( len_trim(crgmask) > 0 ) then
      call atommask( natom, nres, 0, ih(m04), ih(m06), &
         ix(i02), ih(m02), x(lcrd), crgmask, crggp )
      write(6,'(a,a,a,i5,a)') 'Zero-Charge Mask ',crgmask(1:len_trim(crgmask)), ' matches ',sum(crggp(1:natom)),' atoms'
      call remove_charges(crggp, natom, x(l15))
   end if
#endif

   konst = .false.
   belly = .false.

   if( ibelly > 0 .and. (igb > 0 .or. ipb /= 0) ) then

      !          ---here, the only allowable belly has just the first
      !             NATBEL atoms in the moving part.  Check to see that this
      !             requirement is satisfied:

      do i=natbel+1,natom
         if( ix(ibellygp+i-1) /= 0 ) then
            write(6,*) 'When igb>0, the moving part must be at the'
            write(6,*) '   start of the molecule.  This does not seem'
            write(6,*) '   to be the case here.'
            write(6,*) 'natbel,i,igroup(i) = ' &
                  ,natbel,i,ix(ibellygp+i-1)
            call mexit(6,1)
         end if
      end do
   end if


   !     ----- CALCULATE THE SQUARE OF THE BOND PARAMETERS FOR SHAKE
   !           THE PARAMETERS ARE PUT SEQUENTIALLY IN THE ARRAY CONP -----

   do i=1,nbonh + nbona + nbper
      j = ix(iicbh+i-1)
      x(l50+i-1) = req(j)**2
   end do

#ifdef MPI
      if( icfe /= 0 ) then

      !  use the masses of the prmtop file for the first group for both groups:
      !  [only the master nodes communicate here, since non-master nodes
      !   have not yet allocated space]
      ! This leads to problems for dual topology runs, and is therefore skipped
      ! if ifsc is set to one, the masses from both prmtop files are used
         if (ifsc == 0) then
            call mpi_bcast(x(lmass),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
            call mpi_bcast(x(lwinv),natom,MPI_DOUBLE_PRECISION,0,commmaster,ierr)
         end if
         tmass = sum(x(lmass:lmass+natom-1))
         tmassinv = 1.d0/tmass

      !  next, do a minimal sanity check that the SHAKE parameters are
      !  consistent on the two processors:

         ! For Softcore this might be allowed
         ! Put a better check here later
         if( ntc == 2 .and. ifsc == 0) then
            partner = ieor(masterrank,1)
            call mpi_sendrecv( nbonh, 1, MPI_INTEGER, partner, 5, &
                               nbonh_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            call mpi_sendrecv( num_noshake, 1, MPI_INTEGER, partner, 5, &
                               num_noshake_c, 1, MPI_INTEGER, partner, 5, &
                               commmaster, ist, ierr )
            if (qmmm_nml%ifqnt .and. qmmm_nml%qmshake == 0) then
               ! qtw - if qmshake=0, we need to check the QM atoms
               call mpi_sendrecv( &
                       qmmm_struct%nquant, 1, MPI_INTEGER, partner, 5, &
                       nquant_c, 1, MPI_INTEGER, partner, 5, &
                       commmaster, ist, ierr)
               call mpi_sendrecv( &
                       qmmm_struct%noshake_overlap, 1, MPI_INTEGER,partner, 5, &
                       noshake_overlap_c, 1, MPI_INTEGER, partner, 5, &
                       commmaster, ist, ierr)
               if ( (qmmm_struct%nquant-qmmm_struct%noshake_overlap) /= &
                    (nquant_c-noshake_overlap_c) ) then
                  call sander_bomb('mdread2', &
                       'QMMM: NOSHAKE lists are not match in two groups!', &
                       'try noshakemask in cntrl to match the noshake list')
               end if
            else if( nbonh - num_noshake /= nbonh_c - num_noshake_c ) then
               write(6,*) 'SHAKE lists are not compatible in the two groups!'
               call mexit(6,1)
            end if
         else if( ntc == 3 ) then
            write(6,*) 'ntc = 3 is not compatible with icfe>0'
            call mexit(6,1)
         end if

      end if
#endif

      if( igb == 0 .and. ipb == 0 ) &
         call init_extra_pts( &
         ix(iibh),ix(ijbh),ix(iicbh), &
         ix(iiba),ix(ijba),ix(iicba), &
         ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(i40),ix(i42),ix(i44),ix(i46),ix(i48), &
         ix(i50),ix(i52),ix(i54),ix(i56),ix(i58), &
         ih(m06),ix,x,ix(i08),ix(i10), &
         nspm,ix(i70),tmass,tmassinv,x(lmass),x(lwinv),req)

   !  DEBUG input; force checking
#ifdef API
   call load_debug
#else
   call load_debug(5)
#endif

   return
   ! -------------------------------------------------------------------------
   ! Standard format statements:

#ifndef API
   9328 format(/80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)
   9408 format(/4x,'LOADING THE CONSTRAINED ATOMS AS GROUPS',/)
   9409 format(/4x,'LOADING THE TARGETED MD ATOMS AS GROUPS',/)
   9418 format(/4x,'LOADING THE BELLY ATOMS AS GROUPS',/)
   9428 format(/4x,'LOADING THE DECOMP ATOMS AS GROUPS',/)
   9008 format(a80)
#endif

#undef FATAL_ERROR
#undef DELAYED_ERROR
#undef API
end subroutine api_mdread2

! Sets the atomic positions
subroutine set_positions(positions)

   use memory_module, only: natom, lcrd, x

   implicit none

   double precision, intent(in), dimension(natom*3) :: positions

   integer :: i, start

   if (.not. is_setup_) return ! Not sure how to warn...

   start = lcrd - 1
   do i = 1, natom * 3
      x(start + i) = positions(i)
   end do

end subroutine set_positions

! Sets the periodic box vectors
subroutine set_box(a, b, c, alpha, beta, gamma)

    implicit none

    double precision, intent(in) :: a, b, c, alpha, beta, gamma

    if (.not. is_setup_) return ! Not sure how to warn...

    call fill_ucell(a, b, c, alpha, beta, gamma)

end subroutine set_box

subroutine get_box(a_, b_, c_, alpha_, beta_, gamma_)

    use nblist, only : a, b, c, alpha, beta, gamma

    implicit none

    double precision, intent(out) :: a_, b_, c_, alpha_, beta_, gamma_

    if (.not. is_setup_) return ! not sure how to warn

    a_ = a
    b_ = b
    c_ = c
    alpha_ = alpha
    beta_ = beta
    gamma_ = gamma

end subroutine get_box

! Returns the currently active positions
subroutine get_positions(positions)

    use memory_module, only: natom, lcrd, x

    implicit none

    double precision, intent(out), dimension(natom*3) :: positions
    integer :: start, i

    positions(:) = 0.d0
    if(.not. is_setup_) return

    start = lcrd - 1
    do i = 1, natom * 3
        positions(i) = x(start + i)
    end do

end subroutine get_positions

! Computes the energies and forces with the current positions
subroutine energy_forces(energy, forces)

   use memory_module
   use state

   implicit none

   ! Input parameters
   type(potential_energy_rec), intent(out) :: energy
   double precision, dimension(natom*3), intent(out) :: forces

   ! Private variables
   type(state_rec) :: ener
   logical :: qsetup, do_list_update
   integer :: nstep, i, start

   ! Initialization
   energy = null_potential_energy_rec
   forces(:) = 0.d0
   qsetup = .true.
   do_list_update = .true.
   nstep = 1

   if (.not. is_setup_) return ! Not sure how to warn...

   call force(x,ix,ih,ipairs,x(lcrd),x(lforce),ener,ener%vir, &
              x(l96),x(l97),x(l98),x(l99), qsetup, &
              do_list_update,nstep)

   start = lforce - 1
   do i = 1, natom*3
      forces(i) = x(start+i)
   end do
   energy = ener%pot

   return

end subroutine energy_forces

subroutine sander_cleanup()

   use qmmm_adaptive_module, only: adaptive_reset
   use charmm_mod, only : charmm_active, charmm_deallocate_arrays
#ifdef LES
   use genbornles, only: deallocate_gb
#else
   use genborn, only: deallocate_gb
#endif
   use memory_module, only: memory_free
   use molecule, only : deallocate_molecule
   use nblist, only: nblist_deallocate
   use parms, only: clean_parms
   use qmmm_module, only: deallocate_qmmm, qmmm_nml, qmmm_struct, qmmm_vsolv, &
                          qm2_params, qmewald
   use stack, only: deallocate_stacks

   implicit none
#include "../include/md.h"
#include "ew_cntrl.h"

   if (.not. is_setup_) return ! Not sure how to warn...

   call memory_free
   call clean_parms
   if (igb /= 0 .and. igb /= 10 .and. ipb == 0) call deallocate_gb
   call deallocate_stacks
   call nblist_deallocate
   if (allocated(ipairs)) deallocate(ipairs)
   ! periodic does not get reset when ntb==0 inside sander, so do it here
   periodic = 0
   ! QM/MM cleanup
   if (qmmm_nml%ifqnt) then
      call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
      call get_qm2_forces_reset
      qmmm_nml%ifqnt = .false.
      ! Initialize all of the first call flags
      qmewald%ewald_startup = .true.
      qmmm_struct%qm_mm_first_call = .true.
      qmmm_struct%fock_first_call = .true.
      qmmm_struct%fock2_2atm_first_call  = .true.
      qmmm_struct%qm2_allocate_e_repul_first_call = .true.
      qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
      qmmm_struct%qm2_scf_first_call = .true.
      qmmm_struct%zero_link_charges_first_call = .true.
      qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
   end if
   ! CHARMM stuff
   if (charmm_active) call charmm_deallocate_arrays()
   charmm_active = .false.
   ! Reset all of the "first"
   call mdeng_reset
   call mdwrit_reset
   call minwrit_reset
   call pcshift_reset
   call adaptive_reset
   call shake_reset
   call ewald_self_reset
   call nxtsec_reset
   call deallocate_molecule

   is_setup_ = .false.

end subroutine sander_cleanup

! Returns the number of atoms defined in the current system
subroutine sander_natom(n)

   implicit none

#include "../include/memory.h"

   integer, intent(out) :: n

   if (.not. is_setup_) then
      n = 0
   else
      n = natom
   end if

   return

end subroutine sander_natom

! Parses an inpcrd file and fills the coordinate and box arrays with the values
! in the inpcrd file. This routine assumes that the caller already knows how
! many atoms are present and has pre-allocated enough space in the coordinate
! array to hold all of the coordinates. The box variable should have 6 elements.
! If no box is present, then all box lengths and angles will be set to 0. If
! there was any error reading the inpcrd file, ierr will be set to 1.
! Otherwise, it will be 0 on success
subroutine read_inpcrd_file(filename, coordinates, box, ierr)

   use AmberNetcdf_mod
   use binrestart
   use file_io_dat, only : INPCRD_UNIT

   implicit none

   ! Passed variables

   character(len=*), intent(in)  :: filename
   double precision, intent(out) :: coordinates(*)
   double precision, intent(out) :: box(6)
   integer, intent(out)          :: ierr

   ! Local variables

   integer           :: ncid
   integer           :: alloc_failed
   integer           :: natom
   integer           :: i
   character(len=80) :: title
   double precision  :: time
   double precision, dimension(1) :: remd_values
   double precision, allocatable  :: velocities(:)

   box(:) = 0.d0
   ierr = 0

   ! Make sure the NetCDF module doesn't scream at us...
   verbose_netcdf = .false.

   if (NC_checkRestart(filename)) then
      ! See if there are box coordinates
      if (NC_openRead(filename, ncid)) then
         go to 666 ! error label
      end if
      if (NC_readRestartBox(ncid,box(1),box(2),box(3),box(4),box(5),box(6))) &
         box(:) = 0.d0 ! There _is_ no box...
      call NC_close(ncid)
      call read_nc_restart(filename, title, 1, natom, coordinates, velocities, &
                           remd_values, 1, time)
   else
      ! Try reading the raw file
      open(unit=INPCRD_UNIT, file=filename, status='OLD', form='FORMATTED', &
           iostat=alloc_failed)
      if (alloc_failed /= 0) then
         ierr = 1
         return
      end if

      read(INPCRD_UNIT, '(a80)') title
      read(INPCRD_UNIT, *, err=666) natom
      allocate(velocities(natom*3), stat=alloc_failed)
      if (alloc_failed == 1) then
         coordinates(1:natom*3) = 0.d0
         box(:) = 0.d0
         ierr = 1
         return
      end if
      read(INPCRD_UNIT, '(6f12.7)', end=666, err=666) &
            (coordinates(i), i=1,natom*3)
      ! Now get the next 6 numbers. If these are the last numbers in the file,
      ! then we know it is really the set of box dimensions. If there are more
      ! numbers but not enough, it's a corrupt file, otherwise, read all of the
      ! velocities
      read(INPCRD_UNIT, '(6f12.7)', end=667, err=666) &
            (velocities(i), i=1,6)
      read(INPCRD_UNIT, '(6f12.7)', end=668, err=666) &
            (velocities(i), i=7,12)
      read(INPCRD_UNIT, '(6f12.7)', end=666, err=666) &
            (velocities(i), i=13,natom*3)
      ! Now see if we can read our box
      read(INPCRD_UNIT, '(6f12.7)', end=667, err=666) &
            (box(i), i=1,6)
      deallocate(velocities)
      close(INPCRD_UNIT)
   end if

   return

666 continue
   ierr = 1
   close(INPCRD_UNIT)
   if (allocated(velocities)) deallocate(velocities)
   return

667 continue
   ! We come here if we hit EOF trying to read the box, so the box must not
   ! exist. Just deallocate our velocities and leave
   box(:) = 0.d0
   close(INPCRD_UNIT)
   if (allocated(velocities)) deallocate(velocities)
   return

668 continue
   ! We come here if we read in 6 velocities but could not read any more. That
   ! means that the velocities were really the box. Copy them over to the box
   ! and bail out
   close(INPCRD_UNIT)
   if (allocated(velocities)) then
      box(:) = velocities(1:6)
      deallocate(velocities)
   else
      ierr = 1  ! This should never happen
   end if
   return

end subroutine read_inpcrd_file

! This subroutine extracts the number of atoms defined in an inpcrd file. This
! can be used to protect against buffer overruns. natom is set to -1 upon errors
subroutine get_inpcrd_natom(filename, natom)

   use AmberNetcdf_mod
   use binrestart
   use file_io_dat, only : INPCRD_UNIT

   implicit none

   ! Passed variables

   character(len=*), intent(in) :: filename
   integer, intent(out)         :: natom

   ! Local variables

   integer           :: ncid
   integer           :: alloc_failed
   character(len=80) :: title
   double precision  :: time
   integer           :: id1, id2, id3

   ! Make sure the NetCDF module doesn't scream at us...
   verbose_netcdf = .false.

   if (NC_checkRestart(filename)) then
      ! Get the number of atoms
      if (NC_openRead(filename, ncid)) then
         natom = -1
         return
      end if
      if (NC_setupRestart(ncid, title, natom, id1, id2, id3, time)) &
         natom = -1
      call NC_close(ncid)
      return
   else
      ! Try reading the raw file
      open(unit=INPCRD_UNIT, file=filename, status='OLD', form='FORMATTED', &
           iostat=alloc_failed)
      if (alloc_failed /= 0) then
         natom = -1
         return
      end if

      read(INPCRD_UNIT, '(a80)', end=666) title
      read(INPCRD_UNIT, *, err=666, end=666) natom
      close(INPCRD_UNIT)
   end if

   return

666 continue
   natom = -1
   close(INPCRD_UNIT)
   return

end subroutine get_inpcrd_natom

#ifdef LES
end module sanderles_api
#else
end module sander_api
#endif

! Below here we provide the functions _outside_ the module so that the names do
! not become name-mangled (and therefore compiler-dependent). These are just
! wrappers around the module subroutines

! Sets up input options for a gas-phase (i.e., no periodic boundaries) system
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
! gb : integer (optional)
!     If GB model is desired, igb will be set to this.  If not present or 0, a
!     gas phase calculation is set up. If set to 10, a PB calculation is set up
!
! If gb is given, the given GB model is set up. Otherwise, the system is set up
! for a gas-phase simulation

#ifdef LES
#  define SANDER_API_MOD sanderles_api
#else
#  define SANDER_API_MOD sander_api
#endif

subroutine ext_gas_sander_input(inp, gb)

   use SANDER_API_MOD, only : mod_func => gas_sander_input, sander_input

   implicit none

   type(sander_input), intent(out) :: inp
   integer, optional, intent(in)   :: gb

   if (present(gb)) then
      call mod_func(inp, gb)
   else
      call mod_func(inp)
   end if

end subroutine ext_gas_sander_input

! Sets up input options for a periodic system treated with PME
!
! Parameters
! ----------
! inp : type(sander_input)
!     struct of input options that will be filled by this subroutine
subroutine ext_pme_sander_input(inp)

   use SANDER_API_MOD, only : mod_func => pme_sander_input, sander_input

   implicit none

   type(sander_input), intent(out) :: inp

   call mod_func(inp)

end subroutine ext_pme_sander_input

! Initializes a struct with QM options to all default values
!
! Parameters
! ----------
! inp : type(qmmm_input_options)
!     struct of QM input options that will be filled by this subroutine
subroutine ext_qm_sander_input(inp)

   use SANDER_API_MOD, only : mod_func => qm_sander_input, qmmm_input_options

   implicit none

   type(qmmm_input_options), intent(out) :: inp

   call mod_func(inp)

end subroutine ext_qm_sander_input

! Initializes the major data structures needed to evaluate energies and forces
!
! Parameters
! ----------
! prmname : string
!     Name of the topology file
! coordinates : double precision(3*natom)
!     The starting coordinates
! box : double precision(6)
!     The box dimensions
! input_options : type(sander_input)
!     struct of input options used to set up the calculation
! qmmm_options : type(qmmm_input_options) [optional]
!     struct of input options used to set up the QM/MM part of the calculation
! ierr : integer
!     Set to 0 for success or 1 if the setup failed
subroutine ext_sander_setup(prmname, coordinates, box, input_options, qmmm_options, ierr)

   use SANDER_API_MOD, only : mod_func => sander_setup, sander_input, &
                              qmmm_input_options, MAX_FN_LEN

   implicit none

   ! Input parameters
   character(len=MAX_FN_LEN) :: prmname
   double precision, dimension(6), intent(in) :: box
   double precision, dimension(*), intent(in) :: coordinates
   type(sander_input) :: input_options
   type(qmmm_input_options), optional :: qmmm_options
   integer, intent(out) :: ierr

   call mod_func(prmname, coordinates, box, input_options, qmmm_options, ierr)

   return

end subroutine ext_sander_setup

! Sets the atomic positions
subroutine ext_set_positions(positions)

   use memory_module, only: natom
   use SANDER_API_MOD, only : mod_func => set_positions

   implicit none

   double precision, intent(in), dimension(natom*3) :: positions

   call mod_func(positions)

   return

end subroutine ext_set_positions

! Gets the periodic box dimensions
subroutine ext_get_box(a, b, c, alpha, beta, gamma)

    use SANDER_API_MOD, only : mod_func => get_box

    implicit none

    double precision, intent(out) :: a, b, c, alpha, beta, gamma

    call mod_func(a, b, c, alpha, beta, gamma)

    return

end subroutine ext_get_box

! Sets the periodic box vectors
subroutine ext_set_box(a, b, c, alpha, beta, gamma)

   use SANDER_API_MOD, only : mod_func => set_box

   implicit none

   double precision, intent(in) :: a, b, c, alpha, beta, gamma

   call mod_func(a, b, c, alpha, beta, gamma)

   return

end subroutine ext_set_box

! Gets the currently active positions
subroutine ext_get_positions(positions)

    use memory_module, only: natom
    use SANDER_API_MOD, only: mod_func => get_positions

    implicit none

    double precision, intent(out), dimension(natom*3) :: positions

    call mod_func(positions)

    return

end subroutine ext_get_positions

! Computes the energies and forces with the current positions
subroutine ext_energy_forces(energy, forces)

   use SANDER_API_MOD, only : potential_energy_rec, mod_func => energy_forces
   use memory_module, only : natom

   implicit none

   ! Input parameters
   type(potential_energy_rec), intent(out) :: energy
   double precision, dimension(natom*3), intent(out) :: forces

   call mod_func(energy, forces)

   return

end subroutine ext_energy_forces

subroutine ext_sander_cleanup()

   use SANDER_API_MOD, only : mod_func => sander_cleanup

   implicit none

   call mod_func

   return

end subroutine ext_sander_cleanup

! Returns the number of atoms defined in the current system
subroutine ext_sander_natom(n)

    use SANDER_API_MOD, only : mod_func => sander_natom

    implicit none

    integer, intent(out) :: n

    call mod_func(n)

    return

end subroutine ext_sander_natom

! Parses the input coordinate file and properly determines whether velocities
! and/or box are present
subroutine ext_read_inpcrd_file(filename, coordinates, box, ierr)

   use SANDER_API_MOD, only : mod_func => read_inpcrd_file, MAX_FN_LEN

   implicit none

   character(len=MAX_FN_LEN), intent(in) :: filename
   double precision, intent(out) :: coordinates(*)
   double precision, intent(out) :: box(6)
   integer, intent(out)          :: ierr

   call mod_func(filename, coordinates, box, ierr)

   return

end subroutine ext_read_inpcrd_file

! Gets the number of atoms from an inpcrd file
subroutine ext_get_inpcrd_natom(filename, natom)

   use SANDER_API_MOD, only : mod_func => get_inpcrd_natom, MAX_FN_LEN

   implicit none

   character(len=MAX_FN_LEN), intent(in) :: filename
   integer, intent(out) :: natom

   call mod_func(filename, natom)

   return

end subroutine ext_get_inpcrd_natom

! Read a topology file and fill a prmtop_struct with it
subroutine ext_read_prmtop_file(filename, parm, ierr)

    use SANDER_API_MOD, only: mod_func => read_prmtop_file, prmtop_struct, &
                              MAX_FN_LEN

    implicit none

    character(len=MAX_FN_LEN), intent(in) :: filename
    type(prmtop_struct)  :: parm
    integer, intent(out) :: ierr

    ierr = 0

    call mod_func(filename, parm, ierr)

end subroutine ext_read_prmtop_file

! Sets setup to 0 if sander is not set up and 1 if it is
subroutine ext_is_setup(setup)

   use SANDER_API_MOD, only: is_setup

   implicit none

   integer, intent(out) :: setup

   setup = 0
   if (is_setup()) setup = 1

   return

end subroutine ext_is_setup
