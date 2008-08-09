# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 - 2008 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# MolMod is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from molmod.data.periodic import periodic
from molmod.units import A

import os


__all__ = ["mkinput", "mkinput_multiopt"]


template="""%%chk=%(basename)s.chk
%%nproc=%(nproc)i
%%mem=%(mem)s
# %(lot)s %(route_args)s maxdisk=%(maxdisk)s NoSymm

Who cares about the title?

%(charge)i %(spin)i
%(atom_lines)s

%(post)s


"""


def mkinput(
    molecule, charge, spin, lot, route_args, post, nproc, mem, maxdisk, com_filename,
    center=True, overwrite=False, ghost_mask=None
):
    destdir = os.path.dirname(com_filename)
    basename = os.path.basename(com_filename)
    if basename.endswith(".com"):
        basename = basename[:-4]
    if not os.path.isdir(destdir):
        os.makedirs(destdir)
    com_filename = os.path.join(destdir, "%s.com" % basename)
    if not os.path.isfile(com_filename) or overwrite:
        if molecule is None:
            atom_lines = "${atom_lines}"
        else:
            if center:
                # move the coordinates to the origin
                molecule.coordinates -= molecule.coordinates.mean(0)
            symbols = [periodic[number].symbol for number in molecule.numbers]
            # Optionally set ghost atoms:
            if ghost_mask is not None:
                for i in xrange(len(symbols)):
                    if ghost_mask[i]:
                        symbols[i] = "%s-Bq" % symbols[i]
            atom_lines = "\n".join("% 5s % 12.7f % 12.7f % 12.7f" % (
                symbol, cor[0], cor[1], cor[2]
            ) for symbol, cor in zip(symbols, molecule.coordinates/A))
            # Write an xyz file
            molecule.write_to_file(os.path.join(destdir, "geom.xyz"))
        f = file(com_filename, "w")
        f.write(template % {
            "basename": basename,
            "nproc": nproc,
            "mem": mem,
            "lot": lot,
            "route_args": route_args,
            "maxdisk": maxdisk,
            "charge": charge,
            "spin": spin,
            "atom_lines": atom_lines,
            "post": post,
        })
        f.close()



template_multiopt_top="""%%chk=%(basename)s.chk
%%nproc=%(nproc)i
%%mem=%(mem)s
# %(lot)s opt=ModRedundant maxdisk=%(maxdisk)s NoSymm

Who cares about the title?

%(charge)i %(spin)i
%(atom_lines)s

%(post)s


"""

template_multiopt_link="""--Link1--
%%chk=%(basename)s.chk
%%mem=%(mem)s
%%nproc=%(nproc)s
#p %(lot)s opt Geom(AllCheck) maxdisk=%(maxdisk)s NoSymm

Who cares about the title?



"""



def mkinput_multiopt(
    molecule, charge, spin, lot_mem_pairs, post, nproc, maxdisk, com_filename,
    center=True, overwrite=False
):
    destdir = os.path.dirname(com_filename)
    basename = os.path.basename(com_filename)
    if basename.endswith(".com"):
        basename = basename[:-4]
    if len(destdir) > 0 and not os.path.isdir(destdir):
        os.makedirs(destdir)
    com_filename = os.path.join(destdir, "%s.com" % basename)
    if not os.path.isfile(com_filename) or overwrite:
        if center:
            # move the coordinates to the origin
            molecule.coordinates -= molecule.coordinates.mean(0)
        symbols = [periodic[number].symbol for number in molecule.numbers]
        f = file(com_filename, "w")
        # Write a gaussian file (top)
        atom_lines = "\n".join("% 2s % 12.7f % 12.7f % 12.7f" % (
            periodic[number].symbol, cor[0], cor[1], cor[2]
        ) for number, cor in zip(molecule.numbers, molecule.coordinates/A))
        lot, mem = lot_mem_pairs[0]
        f.write(template_multiopt_top % {
            "basename": basename,
            "nproc": nproc,
            "mem": mem,
            "lot": lot,
            "maxdisk": maxdisk,
            "charge": charge,
            "spin": spin,
            "atom_lines": atom_lines,
            "post": post,
        })
        for lot, mem in lot_mem_pairs[1:]:
            f.write(template_multiopt_link % {
                "basename": basename,
                "nproc": nproc,
                "mem": mem,
                "lot": lot,
                "maxdisk": maxdisk,
                "post": post,
            })
        f.close()
        # Write an xyz file
        molecule.write_to_file(os.path.join(destdir, "geom.xyz"))



