#!/bin/sh

#   cdl refinement

if [ "$#" -lt 4 ]; then
   echo "Usage:  phenix.cdl.sh <pdbfile> <mtzfile> <id> <serial> <ciffiles>"
   exit 1
fi

echo "inside phenix.cdl.sh"
echo "$1 $2 $3 $4 $5"

cat <<EOF > ${3}_00$4.eff
refinement {
  input {
    xray_data {
      outliers_rejection = True
      r_free_flags {
        generate = False
      }
    }
  }
  output {
    prefix = "$3"
    serial = $4
    write_eff_file = False
    write_geo_file = False
    write_def_file = False
    write_model_cif_file = False
    write_map_coefficients = True
    export_final_f_model = False
  }
  refine {
    strategy = individual_sites individual_sites_real_space rigid_body \
               individual_adp group_adp tls occupancies group_anomalous
    sites {
    }
    adp { 
      individual {
         anisotropic = none
      }
    }
  }
  main {
    nqh_flips = True
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
    ignore_hydrogens = True
  }
  structure_factors_and_gradients_accuracy {
    algorithm = *fft direct
  }
  gui {
    skip_rsr = True
    skip_kinemage = True
  }
}
EOF

phenix.refine  $1  $2 $5  ${3}_00$4.eff --overwrite

