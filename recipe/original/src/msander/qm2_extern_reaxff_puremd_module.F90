! This is a .F03 file.
!
! module for interfacing with PuReMD code (ReaxFF+EEM in QM/MM mode)
module qm2_extern_reaxff_puremd_module
        use, intrinsic :: iso_c_binding
        implicit none

        private

        interface
                type(c_ptr) function setup_qmmm &
                        (num_qm_atoms, qm_symbols, qm_pos, &
                        num_mm_atoms, mm_symbols, mm_pos_q, sim_box_info, &
                        ffield_filename, control_filename) &
                        bind(C, name='setup_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        integer (c_int), value :: num_qm_atoms
                        type(c_ptr), value :: qm_symbols
                        type(c_ptr), value :: qm_pos
                        integer (c_int), value :: num_mm_atoms
                        type(c_ptr), value :: mm_symbols
                        type(c_ptr), value :: mm_pos_q
                        type(c_ptr), value :: sim_box_info
                        type(c_ptr), value :: ffield_filename
                        type(c_ptr), value :: control_filename
                end function setup_qmmm

                integer (c_int) function reset_qmmm &
                        (handle, num_qm_atoms, qm_symbols, qm_pos, &
                        num_mm_atoms, mm_symbols, mm_pos_q, sim_box_info, &
                        ffield_filename, control_filename) &
                        bind(C, name='reset_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: num_qm_atoms
                        type(c_ptr), value :: qm_symbols
                        type(c_ptr), value :: qm_pos
                        integer (c_int), value :: num_mm_atoms
                        type(c_ptr), value :: mm_symbols
                        type(c_ptr), value :: mm_pos_q
                        type(c_ptr), value :: sim_box_info
                        type(c_ptr), value :: ffield_filename
                        type(c_ptr), value :: control_filename
                end function reset_qmmm

                integer (c_int) function simulate &
                        (handle) &
                        bind(C, name='simulate') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                end function simulate

                integer (c_int) function cleanup &
                        (handle) &
                        bind(C, name='cleanup') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                end function cleanup

                integer (c_int) function set_control_parameter &
                        (handle, keyword, values) &
                        bind(C, name='set_control_parameter') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: keyword
                        type(c_ptr) :: values
                end function set_control_parameter

                integer (c_int) function set_output_enabled &
                        (handle, is_enabled) &
                        bind(C, name='set_output_enabled') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        integer (c_int), value :: is_enabled
                end function set_output_enabled

                integer (c_int) function get_atom_forces_qmmm &
                        (handle, qm_f, mm_f) &
                        bind(C, name='get_atom_forces_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: qm_f
                        type(c_ptr), value :: mm_f
                end function get_atom_forces_qmmm

                integer (c_int) function get_atom_charges_qmmm &
                        (handle, qm_q, mm_q) &
                        bind(C, name='get_atom_charges_qmmm') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: qm_q
                        type(c_ptr), value :: mm_q
                end function get_atom_charges_qmmm

                integer (c_int) function get_system_info &
                        (handle, e_potential, e_kinetic, e_total, temperature, &
                        volume, pressure) &
                        bind(C, name='get_system_info') 
                        use, intrinsic :: iso_c_binding
                        implicit none
                        type(c_ptr), value :: handle
                        type(c_ptr), value :: e_potential
                        type(c_ptr), value :: e_kinetic
                        type(c_ptr), value :: e_total
                        type(c_ptr), value :: temperature
                        type(c_ptr), value :: volume
                        type(c_ptr), value :: pressure
                end function get_system_info
        end interface

        public :: get_reaxff_puremd_forces, reaxff_puremd_finalize

        character(len=*), parameter, public :: module_name = "qm2_extern_reaxff_puremd_module"
        type(c_ptr), save, private :: handle = C_NULL_PTR
        
        type reaxff_nml_type
             character(len=20) :: control        ! control file name
             character(len=20) :: ffield         ! force field file name
             integer           :: numthreads    ! number of threads for reaxff-puremd (OpenMP only)
             real              :: solvtol        ! charge solver tolerance
             integer           :: solvmaxit      ! max number of iterations for the charge solver
             integer           :: solvprecond    ! preconditioner type for the charge solver
             real              :: nbrcut         ! bonded interaction cutoff
             real              :: hbondcut       ! hydrogen bonding interaction cutoff
             real              :: thbcut         ! valence cutoff (based on bond-order)
        end type reaxff_nml_type


contains

        subroutine get_reaxff_puremd_forces( num_qm_atoms, qm_pos, qm_types, &
                        qm_q, num_mm_atoms, mm_pos_q, mm_types, e_total, &
                        qm_f, mm_f, qmcharge )
                use, intrinsic :: iso_c_binding
                use ElementOrbitalIndex, only : elementSymbol
                implicit none

                integer (c_int), intent(in) :: num_qm_atoms                     ! number of QM atoms
                real (c_double), target, intent(in) :: qm_pos(3,num_qm_atoms)   ! QM atom coordinates
                integer (c_int), target, intent(in) :: qm_types(num_qm_atoms)   ! QM atom types
                real (c_double), target, intent(inout) :: qm_q(num_qm_atoms)    ! QM atom charges (nuclear charge in au)
                integer (c_int), intent(in) :: num_mm_atoms                     ! number of MM atoms
                real (c_double), target, intent(in) :: mm_pos_q(4,num_mm_atoms) ! MM atom coordinates and charges (nuclear charge in au)
                integer (c_int), target, intent(in) :: mm_types(num_mm_atoms)   ! MM atom types
                real (c_double), target, intent(out) :: e_total                 ! SCF energy (kcal/mol)
                real (c_double), target, intent(out) :: qm_f(3,num_qm_atoms)    ! SCF QM force (AMU * Angstroms / ps^2)
                real (c_double), target, intent(out) :: mm_f(3,num_mm_atoms)    ! SCF MM force (AMU * Angstroms / ps^2)
                integer (c_int), intent(in) :: qmcharge                         ! total charge of the QM region

                logical, save :: first_call = .true.
                type(reaxff_nml_type), save :: reaxff_nml
                logical :: control_exists
                logical :: ffield_exists
                integer :: ix
                integer (c_int) :: ret
                character(kind=c_char, len=1024), target :: control_filename, ffield_filename, keyword, values
                ! triplets for lengths and angles of QM region simulation box (Angstroms and degrees)
                real (c_double), target :: sim_box_info(6)
                character(len=2), dimension(num_qm_atoms), target :: qm_symbols
                character(len=2), dimension(num_mm_atoms), target :: mm_symbols

#if defined(HAVE_REAXFF_PUREMD)
                ! NOTE: PuReMD must run with periodic boundary conditions (PBCs) ON,
                !       so to compensate the simulation box will have void space added around it
                !       (20 angstroms, as the long-range cut-off is 10 angstroms) in order
                !       negate the effect of PBCs for QMMM
                do ix = 1, 3
                        sim_box_info(ix) = MAX(MAXVAL(mm_pos_q(ix,:)), MAXVAL(qm_pos(ix,:))) &
                                - MIN(MINVAL(mm_pos_q(ix,:)), MINVAL(qm_pos(ix,:))) + 20.0
                end do
                ! orthogonal simulation box
                sim_box_info(4:6) = 90.0

                ! get character representations of atoms
                do ix = 1, num_qm_atoms
                        qm_symbols(ix) = elementSymbol(qm_types(ix))
                end do
                do ix = 1, num_mm_atoms
                        mm_symbols(ix) = elementSymbol(mm_types(ix))
                end do
                if ( first_call ) then
                        first_call = .false.

                        write (6, '(/,a,/)') '  >>> Running calculations with ReaxFF <<<'

                        call get_namelist_reaxff(reaxff_nml)
                        call print_namelist_reaxff(reaxff_nml)

                        control_filename = trim(adjustl(reaxff_nml%control)) // C_NULL_CHAR
                        ffield_filename = trim(adjustl(reaxff_nml%ffield))  // C_NULL_CHAR
 
                        inquire( file=control_filename, exist=control_exists )
                        inquire( file=ffield_filename, exist=ffield_exists )

                        ! check if the force field file exists
                        if ( ffield_exists .eqv. .false. ) then
                                call sander_bomb("get_reaxff_puremd_forces", &
                                        "ERROR: Specified force field file does not exist!", &
                                        "Will exit now")
                        endif

                        ! if reaxff-puremd control supplied by user file
                        ! and exists, use it to set simulation parameters
                        if ( len_trim(reaxff_nml%control) > 0 .and. control_exists ) then
                                handle = setup_qmmm( num_qm_atoms, c_loc(qm_symbols), &
                                        c_loc(qm_pos), num_mm_atoms, c_loc(mm_symbols), &
                                        c_loc(mm_pos_q), c_loc(sim_box_info), &
                                        c_loc(ffield_filename), c_loc(control_filename) )
                        else
                                if ( len_trim(reaxff_nml%control) > 0 ) then
                                    write(6,*) "[WARNING] Speficied control file for ReaxFF does not exist!"
                                endif

                                handle = setup_qmmm( num_qm_atoms, c_loc(qm_symbols), &
                                        c_loc(qm_pos), num_mm_atoms, c_loc(mm_symbols), &
                                        c_loc(mm_pos_q), c_loc(sim_box_info), &
                                        c_loc(ffield_filename), C_NULL_PTR )
                        endif

                        ! NVE ensemble
                        keyword = "ensemble_type" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! MD steps
                        keyword = "nsteps" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! time step length (in fs)
                        keyword = "dt" // C_NULL_CHAR
                        values = "0.25" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! OpenMP number of threads
                        keyword = "num_threads" // C_NULL_CHAR
                        write (values, *) reaxff_nml%numthreads
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! enable periodic boundary conditions
                        keyword = "periodic_boundaries" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! do not remap atom coordinates within simulation box boundaries
                        keyword = "reposition_atoms" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! recompute Verlet neighbor lists at every (1) MD step
                        keyword = "reneighbor" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! disable force and energy tabulation for Coulomb interactions
                        keyword = "tabulate_long_range" // C_NULL_CHAR
                        values = "0" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! calculate energies at every (1) MD step
                        keyword = "energy_update_freq" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! add a 2.5 Angstrom buffer to Verlet neighbor list cut-off
                        keyword = "vlist_buffer" // C_NULL_CHAR
                        values = "2.5" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 5.0 Angstrom bond interaction cut-off
                        keyword = "nbrhood_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%nbrcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 0.001 threshold for valence angle interactions
                        keyword = "thb_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%thbcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 7.5 Angstrom hydrogen bond interaction cut-off
                        keyword = "hbond_cutoff" // C_NULL_CHAR
                        write (values, *) reaxff_nml%hbondcut
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! 0.3 Angstrom bond graph calculation cut-off
                        keyword = "bond_graph_cutoff" // C_NULL_CHAR
                        values = "0.3" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! EEM model (full system) for charge calculations
                        keyword = "charge_method" // C_NULL_CHAR
                        values = "1" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! net charge for system (in Coulombs)
                        keyword = "cm_q_net" // C_NULL_CHAR
                        write (values, *) qmcharge
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! conjugate gradient algorithm in charge solver
                        keyword = "cm_solver_type" // C_NULL_CHAR
                        values = "2" // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! max. iterations in charge solver
                        keyword = "cm_solver_max_iters" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvmaxit
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! tolerance in charge solver
                        keyword = "cm_solver_q_err" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvtol
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! Jacobi preconditioner in charge solver
                        keyword = "cm_solver_pre_comp_type" // C_NULL_CHAR
                        write (values, *) reaxff_nml%solvprecond
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ! disable file I/O
                        ret = set_output_enabled( handle, 0 )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_output_enabled", &
                                "Will exit now")

                        ret = simulate( handle )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::simulate", &
                                "Will exit now")

                        ret = get_atom_forces_qmmm( handle, c_loc(qm_f), c_loc(mm_f) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_forces_qmmm", &
                                "Will exit now")

                        ! disregard MM atom charges, as static (input-only)
                        ret = get_atom_charges_qmmm( handle, c_loc(qm_q), C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_charges_qmmm", &
                                "Will exit now")

                        ! disregard all values except total energy
                        ret = get_system_info( handle, C_NULL_PTR, C_NULL_PTR, &
                                c_loc(e_total), C_NULL_PTR, C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_system_info", &
                                "Will exit now")
                else
                        ret = reset_qmmm( handle, num_qm_atoms, c_loc(qm_symbols), &
                                c_loc(qm_pos), num_mm_atoms, c_loc(mm_symbols), &
                                c_loc(mm_pos_q), c_loc(sim_box_info), &
                                C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::reset_qmmm", &
                                "Will exit now")

                        ! net charge for system (in Coulombs)
                        keyword = "cm_q_net" // C_NULL_CHAR
                        write (values, *) qmcharge
                        values = trim(adjustl(values)) // C_NULL_CHAR
                        ret = set_control_parameter( handle, c_loc(keyword), c_loc(values) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::set_control_parameter", &
                                "Will exit now")

                        ret = simulate( handle )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::simulate", &
                                "Will exit now")

                        ret = get_atom_forces_qmmm( handle, c_loc(qm_f), c_loc(mm_f) )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_forces_qmmm", &
                                "Will exit now")

                        ! disregard MM atom charges, as static (input-only)
                        ret = get_atom_charges_qmmm( handle, c_loc(qm_q), C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_atom_charges_qmmm", &
                                "Will exit now")

                        ! disregard all values except total energy
                        ret = get_system_info( handle, C_NULL_PTR, C_NULL_PTR, &
                                c_loc(e_total), C_NULL_PTR, C_NULL_PTR, C_NULL_PTR )
                        if ( ret /= 0_c_int ) call sander_bomb("get_reaxff_puremd_forces", &
                                "ERROR: get_reaxff_puremd_forces::get_system_info", &
                                "Will exit now")
                end if
#else
                call sander_bomb('get_reaxff_puremd_forces','reaxff-puremd is not enabled', &
                        'Check your installation or reconfigure with the -reaxff-puremd option.')
#endif
        end subroutine get_reaxff_puremd_forces


        subroutine reaxff_puremd_finalize( )
                use, intrinsic :: iso_c_binding
                implicit none

                integer (c_int) :: ret

#if defined(HAVE_REAXFF_PUREMD)
                if ( c_associated(handle) ) then
                        ret = cleanup( handle )
                        if ( ret /= 0_c_int ) call sander_bomb("reaxff_puremd_finalize", &
                                "ERROR: reaxff_puremd_finalize::cleanup", &
                                "Will exit now")
                endif
#else
                call sander_bomb('reaxff_puremd_finalize','reaxff-puremd is not enabled', &
                        'Check your installation or reconfigure with the -reaxff-puremd option.')
#endif
        end subroutine reaxff_puremd_finalize


        subroutine get_namelist_reaxff( reaxff_nml )
                implicit none

                type(reaxff_nml_type), intent(out) :: reaxff_nml

                character(len=20) :: control        ! control file name
                character(len=20) :: ffield         ! force field file name
                integer           :: numthreads     ! number of threads for reaxff-puremd (OpenMP only)
                real              :: solvtol        ! charge solver tolerance
                integer           :: solvmaxit      ! max number of iterations for the charge solver
                integer           :: solvprecond    ! preconditioner type for the charge solver
                real              :: nbrcut         ! bonded interaction cutoff
                real              :: hbondcut       ! hydrogen bonding interaction cutoff
                real              :: thbcut         ! valence cutoff (based on bond-order)
                integer :: ierr

                namelist /reaxff/ control, ffield, numthreads, solvtol, &
                        solvmaxit, solvprecond, nbrcut, hbondcut, thbcut

                ! default values
                control = '' 
                ffield = 'ffield.reaxff'
                numthreads = 1
                solvtol = 1.0e-8
                solvmaxit = 200
                solvprecond = 1 ! Jacobi preconditioner for charge solver
                nbrcut = 5.0 ! Angstroms
                hbondcut = 7.5 ! Angstroms
                thbcut = 0.001 

                ! Read namelist
                rewind 5
                read(5, nml=reaxff, iostat=ierr)
                if ( ierr > 0 ) then
                        call sander_bomb('get_namelist_reaxff (qm2_extern_reaxff_puremd_module)', &
                        '&reaxff namelist read error', &
                        'Please check your input.')
                else if ( ierr < 0 ) then
                        write(6,'(a)') '&reaxff namelist read success'
                end if

                reaxff_nml%control = control
                reaxff_nml%ffield = ffield
                reaxff_nml%numthreads = numthreads
                reaxff_nml%solvtol = solvtol
                reaxff_nml%solvmaxit = solvmaxit
                reaxff_nml%solvprecond = solvprecond
                reaxff_nml%nbrcut = nbrcut
                reaxff_nml%hbondcut = hbondcut
                reaxff_nml%thbcut = thbcut
        end subroutine get_namelist_reaxff


        subroutine print_namelist_reaxff( reaxff_nml )
                implicit none

                type(reaxff_nml_type), intent(in) :: reaxff_nml        

                write(6,'(a)') '---------ReaxFF options-------'
                ! Print the control file name if it is set (advanced option)
                if ( len_trim(reaxff_nml%control) > 0 ) then
                        write(6,'(2a)') ' control      ', reaxff_nml%control
                endif
                write(6,'(2a)')         ' ffield       ', reaxff_nml%ffield
                write(6,'(2a)')         ' solvtype     ', "EEM(fixed)"
                write(6,"(A,E7.1)")     ' solvtol      ', reaxff_nml%solvtol
                write(6,"(A,I7)")       ' solvmaxit    ', reaxff_nml%solvmaxit
                write(6,"(A,I7)")       ' solvprecond  ', reaxff_nml%solvprecond
                write(6,"(A,F7.2)")     ' nbrcut       ', reaxff_nml%nbrcut
                write(6,"(A,F7.2)")     ' hbondcut     ', reaxff_nml%hbondcut
                write(6,"(A,F7.5)")     ' thbcut       ', reaxff_nml%thbcut
                write(6,"(A,I7,A)")     ' numthreads   ', reaxff_nml%numthreads, '(OpenMP only)'
                write(6,'(a)')          '-----end ReaxFF options-------'
        end subroutine print_namelist_reaxff
end module qm2_extern_reaxff_puremd_module
