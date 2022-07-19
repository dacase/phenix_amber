#!/bin/sh

diff amber_00$1.log amber$1.log | awk 'NF==8 && $2!="Amber" {print $2}' \
     > amber_energies$1.dat
