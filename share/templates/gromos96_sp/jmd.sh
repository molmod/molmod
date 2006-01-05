#!/bin/sh
#$Id: jmd.sh $

rm -f fort.*
#---input units
ln -s     mta1.dat         fort.20    # molecular topology
ln -s     sxmd1.dat        fort.21    # initial coordinates

#----output units
ln -s     sxmd2.dat        fort.11    # final coordinates
ln -s     rxmd2.dat        fort.12    # coordinate trajectory

#---run the program
${gromos_root}/promd.64 < imd.dat > omd.out

#---clean up after us
rm -f fort.*
