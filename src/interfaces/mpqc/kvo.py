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

from keyval import KeyValObject

from pychem.moldata import periodic


def create_molecule(molecule):
    return KeyValObject(
        'Molecule',
        {
            'atoms': [periodic.symbol[number] for number in molecule.numbers],
            'geometry': molecule.coordinates.tolist()
        }
    )

  
def create_std_functional(functional):
    return KeyValObject(
        'StdDenFunctional', 
        {'name': functional}
    )


def create_mole(molecule_kvo, charge, method, basis, functional=None, guess=True):
    result = KeyValObject(method, {
        'molecule': molecule_kvo,
        'total_charge': charge,
        'basis': KeyValObject('GaussianBasisSet', {
            'molecule': molecule_kvo,
            'name': basis
        })
    })
    if functional != None:
        result.functional = create_std_functional(functional)
    if guess:
        result.guess_wavefunction = create_mole(molecule_kvo, charge, "CLHF", basis="STO-3G", guess=False)
    return result


def create_single_point(molecule, charge, method, basis, functional=None, guess=True):        
        molecule_kvo = create_molecule(molecule)
        
        return KeyValObject(attributes={
            'molecule': molecule_kvo,
            'mpqc': KeyValObject(attributes={
                'mole': create_mole(molecule_kvo, charge, method, basis, functional, guess)
            })
        })


def create_optimize(molecule, charge, method, basis, functional=None, guess=True):        
        molecule_kvo = create_molecule(molecule)
        mole_kvo = create_mole(molecule_kvo, charge, method, basis, functional, guess)
        
        return KeyValObject(attributes={
            'molecule': molecule_kvo,
            'mpqc': KeyValObject(attributes={
                'optimize': 1,
                'mole': mole_kvo,
                'opt': KeyValObject('QNewtonOpt', {
                    'function': mole_kvo,
                    'update': KeyValObject('BFGSUpdate'),
                    'convergence': KeyValObject('MolEnergyConvergence', {
                        'cartesian': 1,
                        'energy': mole_kvo
                    })
                })
            })
        })
