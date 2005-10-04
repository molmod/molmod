# PyChem is a general chemistry oriented python package.
# Copyright (C) 2005 Toon Verstraelen
# 
# This file is part of PyChem.
# 
# PyChem is free software; you can redistribute it and/or
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

from keyval import KeyVal, KeyValObject

from pychem.moldata import periodic

def create_molecule_kvo(molecule):
    return KeyValObject(
        'Molecule',
        named_items=[
            ('atoms', [periodic.symbol[number] for number in molecule.numbers]),
            ('geometry', molecule.coordinates.tolist())
        ]
    )
    
def create_std_functional_kvo(functional):
    if functional != None:
        return KeyValObject(
            'StdDenFunctional', 
            named_items=[('name', functional)]
        )

def create_guess_kvo(molecule_kvo, charge, method="CLHF"):
    return KeyValObject(method, named_items=[
        ('molecule', molecule_kvo),
        ('total_charge', charge),
        ('basis', KeyValObject("GaussianBasisSet", named_items=[
            ('molecule', molecule_kvo),
            ('name', 'STO-3G')
        ]))
    ])
    
def create_mole_kvo(molecule_kvo, charge, method, basis, functional=''):
    return KeyValObject(method, named_items=[
        ('molecule', molecule_kvo),
        ('total_charge', charge),
        ('basis', KeyValObject('GaussianBasisSet', named_items=[
            ('molecule', molecule_kvo),
            ('name', basis)
        ])),
        ('functional', create_std_functional_kvo(functional)),
        ('guess_wavefunction', create_guess_kvo(molecule_kvo, charge))
    ])

def create_single_point_kv(molecule, charge, method, basis, functional=''):        
        molecule_kvo = create_molecule_kvo(molecule)
        
        return KeyVal(named_items=[
            ('molecule', molecule_kvo),
            ('mpqc', KeyValObject(named_items=[
                ('mole', create_mole_kvo(molecule_kvo, charge, method, basis, functional))
            ]))
        ])       


def create_optimize_kv(molecule, charge, method, basis, functional=''):        
        molecule_kvo = create_molecule_kvo(molecule)
        mole_kvo = create_mole_kvo(molecule_kvo, charge, method, basis, functional)
        
        return KeyVal(named_items=[
            ('molecule', molecule_kvo),
            ('mpqc', KeyValObject(named_items=[
                ('optimize', 1),
                ('mole', mole_kvo),
                ('opt', KeyValObject('QNewtonOpt', named_items=[
                    ('function', mole_kvo),
                    ('update', KeyValObject('BFGSUpdate')),
                    ('convergence', KeyValObject('MolEnergyConvergence', named_items=[
                        ('cartesian', 1),
                        ('energy', mole_kvo)
                    ]))
                ]))
            ]))
        ])
