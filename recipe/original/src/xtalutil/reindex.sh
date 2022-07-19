phenix.mtz.dump -f s -c ../XrayPrep/4yuo.mtz | tr ',' '\t' | \
    awk 'NR>1 {print $1*2,$2*2,$3,$4,$5,$6}' > foo.fmtz

f2mtz hklin foo.fmtz hklout 4yuo.mtz <<EOF
CELL 85.800  104.860   89.110  90.00  90.00  90.00
SYMMETRY P1
LABOUT H K L FOBS SIGFOBS R-free-flags
CTYPOUT H H H F Q I
TITLE 4yuo, reindexed to 2x2x1
END
EOF

/bin/rm -f foo.fmtz
