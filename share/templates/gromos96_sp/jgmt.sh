#!/bin/sh
#$Id: jgmt.sh  $

# Start by removing old fortran links
rm -f fort.*

ln -s     mta1.dat          fort.10    # generated molecular topology
ln -s     mtbZEO.dat        fort.11    # residue topology building blocks
ln -s     ifpZEO.dat        fort.12    # interaction function parameters

${gromos_root}/progmt.64 < igmt.dat > ogmt.out

#clean up
rm -f  fort.*

