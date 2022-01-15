#define STATEINF_FLD_C 5
#define TITR_RES_C 50
#define TITR_STATES_C 200
#define ATOM_CHRG_C 1000
#define MAX_H_COUNT 4

type :: const_e_info
   sequence
   integer :: num_states, first_atom, num_atoms, first_state, first_charge
end type const_e_info

type (const_e_info), parameter :: NULL_CE_INFO = const_e_info(0,0,0,0,0)
