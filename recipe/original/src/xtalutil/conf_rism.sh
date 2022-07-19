#!/bin/sh

# take a single conformer from a solute pdb file ($1, default is dry.pdb)
#     run RISM to get a ccp4 map file of the solvent, suitable for
#     providing to conf_solvent.sh

# generally will require minor editing below


# ================== convert CYS and HIS residues; add symmetry cards
cat s19.symm $1.pdb | sed -e 's/CYS/CYX/' | sed -e 's/HIS/HID/' > dry1.pdb

# ================== get unit cell pdb
UnitCell -p dry1.pdb -o uc.pdb

# ================== set up in Amber
cat <<EOF > leap.in
source oldff/leaprc.ff14SB
set default nocenter on
x = loadpdb uc.pdb
set x box {45.90 40.70 30.10}
saveamberparm x $1.parm7 $1.rst7
savepdb x $1.amber.pdb
quit
EOF

tleap -f leap.in

/bin/rm -f dry1.pdb leap.in

# ============= Run rism3d.snglpnt
solvent=/home/case/projects/3drism/solvent/tienhung/cSPCE/KH/water/rism.xvv

/home/case/bin/rism3d.snglpnt --periodic pme --pdb $1.amber.pdb \
   --prmtop $1.parm7 --rst $1.rst7 \
    --xvv $solvent \
    --closure kh --tolerance 0.000001 --solvcut 9. \
    --grdspc 0.35,0.35,0.35 --buffer 1.0 --verbose 2 --volfmt ccp4 \
    --electronMap $1.khe --centering 0 \
     > $1.kh.r3d

ln -s $1.khe.O.1.ccp4 a0.map

# convert to proper axis order for SFALL
mapmask mapin a0.map mapout $1.map << EOF
AXIS Z X Y
EOF

# CLEAN:
/bin/rm -f a0.map 
