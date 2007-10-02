# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2007 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

import numpy

from scipy.special import gamma


def gaussian(r, p, a):
    # not normalized
    return sum(r**p)*numpy.exp(-a*numpy.dot(r, r))

def normalization(p, a):
    return (
        gamma(0.5*(p[0]+1))*
        gamma(0.5*(p[1]+1))*
        gamma(0.5*(p[2]+1))*
        a**(-1.5-0.5*sum(p))*(
            1 + (-1)**p[0] + (-1)**p[1] + (-1)**p[2] +
            (-1)**(p[0]+p[1]) + (-1)**(p[1]+p[2]) + (-1)**(p[2]+p[0]) +
            (-1)**(p[0]+p[1]+p[2])
    ))/8.0


