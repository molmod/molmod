# MolMod is a collection of molecular modelling tools for python.
# Copyright (C) 2005 Toon Verstraelen
#
# This file is part of MolMod.
#
# MolMod is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# --
"""
A set of tools to calculate tunneling corrections.

This is an example script (or input file as you may want to call it), that
shows how the tunneling module works:

>>> from molmod.units import kcalmol, kjmol, invcm
>>> from molmod.tunneling import *
>>>
>>> factor, error = eckart(T=400, Ef=30.33*kjmol, Er=69.13*kjmol, nu=845.9*invcm)
>>> print factor, error
1.51895388399 8.63859902856e-07
>>> print eckart(T=400, Ef=30.33*kjmol, Er=69.13*kjmol, nu=845.9*invcm)
(1.518953883993303, 8.6385990285569791e-07)
>>> print wigner(T=400, nu=845.9*invcm)
(1.3857382192095034, 0.0)
>>>
>>> import numpy
>>> temperatures = numpy.arange(25, 1101, 250)
>>> print_batch("1-2(H)", temperatures, eckart, Ef=30.33*kjmol, Er=69.13*kjmol, nu=845.9*invcm)
   label    T[K]   eckart-correction       error
    1-2(H)    25    6.594304538e+50    7.494871739e+42
    1-2(H)   275    2.506815637e+00    6.945792866e-06
    1-2(H)   525    1.275367654e+00    1.046503151e-06
    1-2(H)   775    1.121730010e+00    8.536537389e-07
    1-2(H)  1025    1.070370556e+00    6.519441615e-07
>>> print_batch("1-2(H)", temperatures, eckart, Ef=30.33*kjmol, Er=69.13*kjmol, nu=845.9*invcm)
   label    T[K]   eckart-correction       error
    1-2(H)    25    6.594304538e+50    7.494871739e+42
    1-2(H)   275    2.506815637e+00    6.945792866e-06
    1-2(H)   525    1.275367654e+00    1.046503151e-06
    1-2(H)   775    1.121730010e+00    8.536537389e-07
    1-2(H)  1025    1.070370556e+00    6.519441615e-07
>>> print_batch("1-2(H)", temperatures, wigner, nu=845.9*invcm)
   label    T[K]   wigner-correction       error
    1-2(H)    25    9.974898412e+01    0.000000000e+00
    1-2(H)   275    1.816107307e+00    0.000000000e+00
    1-2(H)   525    1.223920599e+00    0.000000000e+00
    1-2(H)   775    1.102756487e+00    0.000000000e+00
    1-2(H)  1025    1.058744190e+00    0.000000000e+00
"""


from molmod.units import unified, kjmol
from molmod.constants import boltzman

from scipy.integrate import quad, Inf
import numpy

__all__ = ["eckart", "wigner", "print_batch"]


h = 2*numpy.pi

def eckart(T, Ef, Er, nu):
    """Computes the Eckart correction factor for the given parameters.

    Arguments
        T  --  the temperature
        Ef  --  the forward energy barrier
        Er  --  the reverse energy barrier
        nu  --  the imaginary frequency (as a real number)

    Returns: factor, error
        factor  --  the Eckart correction factor
        error  -- the error on the correction factor
    """

    l = (Ef**(-0.5) + Er**(-0.5))**(-1)*numpy.sqrt(2) / nu

    def alpha(E):
        return numpy.sqrt(2*l**2*E/h**2)

    def beta(E):
        return numpy.sqrt(2*l**2*(E -(Ef-Er))/h**2)

    def delta(E):
        return numpy.sqrt(4*Ef*Er/(h*nu)**2-0.25)

    def P(E):
        return (
            numpy.cosh(2*numpy.pi*(alpha(E) + beta(E))) -
            numpy.cosh(2*numpy.pi*(alpha(E) - beta(E)))
        ) / (
            numpy.cosh(2*numpy.pi*(alpha(E) + beta(E))) +
            numpy.cosh(2*numpy.pi*(delta(E)))
        )

    def integrandum(E):
        return P(E)*numpy.exp(-(E-Ef)/(boltzman*T))
    emin=max([0, Ef-Er])
    emax=500*kjmol
    energies = numpy.arange(emin, emax, 1*kjmol)
    integranda = numpy.array([integrandum(energy) for energy in energies])
    if max(integranda) * 1e-5 < max([integranda[0], integranda[-1]]):
        print "Integrandum is not negligible at borders.", integranda[0] / max(integranda), integranda[-1] / max(integranda)

    integral, error = quad(integrandum, emin, emax)
    factor = 1.0/(boltzman*T)
    return integral*factor, error*factor
eckart.label = "eckart"


def wigner(T, nu):
    """Returns the Wigner tunneling correction.

    Arguments
        T  --  The temperature
        nu  --  The imaginary frequency

    Returns: factor, error
        factor  --  The Wigner correction factor
        error  --  Always zero
    """

    return 1 + (h*nu/(boltzman*T))**2/24, 0.0
wigner.label = "wigner"


def print_batch(label, temperatures, correction, *args, **kwargs):
    """ A tool that computes and prints correction factors for a temperature range.

    Arguments
        label  --  label that is printed in the first column
        temperatures  --  sequence of temperatures (list, array, ...)
        correction  --  function that computes the correction factor, i.e. wigner
                        or eckart

    After these three arguments, a series of arguments can be given, which will
    be passed to the function that computes the correction factor. (see examples)
    """
    print "   label    T[K] % 8s-correction       error  " % correction.label
    for temperature in temperatures:
        factor, error = correction(temperature, *args, **kwargs)
        print "% 10s % 5.0f   % 12.9e   % 12.9e" % (label, temperature, factor, error)

