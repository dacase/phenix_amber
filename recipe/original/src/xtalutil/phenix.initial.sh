#!/bin/sh

#   analyze input structure

if [ "$#" -ne 3 ]; then
   echo "Usage:  phenix.initial.sh <pdbfile> <mtzfile> <cif-files>"
   exit 1
fi

cat <<EOF > refine.eff
refinement {
  input {
    xray_data {
      outliers_rejection = False
      r_free_flags {
        generate = True
      }
    }
  }
  output {
    prefix = "initial"
    serial = 1
    write_eff_file = True
    write_geo_file = False
    write_def_file = False
    write_model_cif_file = False
    write_map_coefficients = False
    export_final_f_model = True
  }
  refine {
    strategy = individual_sites individual_sites_real_space rigid_body \
               individual_adp group_adp tls occupancies group_anomalous
  }
  main {
    nqh_flips = False
    number_of_macro_cycles = 1
    target = auto *ml mlhl ml_sad ls mli
    use_experimental_phases = False
    scattering_table = wk1995 *it1992 n_gaussian electron neutron
  }
  hydrogens {
    refine = individual *riding Auto
  }
  pdb_interpretation {
    c_beta_restraints = False
  }
  mask {
    ignore_hydrogens = False
  }
  structure_factors_and_gradients_accuracy {
    algorithm = fft *direct
  }
  gui {
    skip_rsr = True
    skip_kinemage = True
  }
}
EOF

phenix.refine  $1  $2 $3  refine.eff --overwrite

/bin/rm -f refine.eff
