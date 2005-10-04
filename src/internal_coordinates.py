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

from pychem.molecular_graphs import MolecularGraph

import math, copy
import Numeric, LinearAlgebra

__all__ = [
    "InternalCoordinate", "Select", 
    "Binary", "Add", "Sub", "Delta", "Mul", "Div", "Dot", "Cos", 
    "Unary", "Distance", "DistanceSqr", "Sqr", "Sqrt", "ArcCos",
    "Configuration", "JacobianSolver"
]


# Classes for calculating the partial derivates of the internal
# coordinates towards the carthesian coordinates

class InternalCoordinate(object):
    """
    This is the base class for the internal coordinates. The following
    attributs are availlable for all internal coordinates. They are all used
    for plotting purposes:
        - name: a (latex) description for the internal coordinate
        - unit: a (latex) description of the unit used for this internal
        coordinate
        - conversion: a conversion factor in case you like to convert this
        variable to another unit
        - ylim: the range for the partial derivate of the molecular energy
        towards this internal coordinate (after conversion of ordinate unit)
        
    An instance of the InternalCoordinate class acts as a callable. It takes
    only one parameter:
        - coordinates: a set of molecular carthesian coordinates
    and returns a tuple with:
        - values: the values of the internal coordinate
        - gradient: the partial derivates of the internal coordinates towards
        the carthesian coordinates. This matrix has the same shape as
        coordinates
        
    For the return values there are two conventions:
        A) The returned value of the internal coordinate is one-dimensional.
            This is the regular case. the gradient contains the regular
            partial derivates.       
        B) The returned value of the internal coordinate is three-dimensional.
            The i'th element of the returned internal coordinates is only 
            dependant on the i'th column of the given coordinates and the
            i'th column of the gradient contains the partial derivates of the
            i'th internal coordinate towards the given carthesian coordinates.
    Some Internal coordinate classes can handle both conventions (eg. Delta)
    """
    
    def __init__(self, **keyvals):
        self.__dict__.update(keyvals)
        
    def symbol(self):
        temp = str(self.__class__)
        return temp[temp.rfind(".")+1:-2].lower()
        
    def label(self):
        raise NotImplementedError


class Select(InternalCoordinate):
    """
    This internal coordinate is only a 'helper' for calculating other
    internal coordinates. It selects a given atom from the molecule.
    """

    def __init__(self, index, **keyvals):
        InternalCoordinate.__init__(self, **keyvals)
        self.index = index
        
    def __call__(self, coordinates):
        gradient = Numeric.zeros(coordinates.shape, Numeric.Float)
        gradient[self.index,:] = 1
        return coordinates[self.index], gradient
        
    def label(self):
        return "%s(%i)" % (self.symbol(), self.index)


class Binary(InternalCoordinate):
    """
    The base class for internal coordinates that take _two_ internal
    coordinate as input parameter. 'iic' in the code stands for "input 
    internal coordinate".
    """
       
    def __init__(self, iic1, iic2, **keyvals):
        InternalCoordinate.__init__(self, **keyvals)
        self.iic1 = iic1
        self.iic2 = iic2
        
    def label(self):
        return "%s(%s,%s)" % (self.symbol(), self.iic1.label(), self.iic2.label())

class Add(Binary):
    def __call__(self, coordinates):
        v1, gv1 = self.iic1(coordinates)
        v2, gv2 = self.iic2(coordinates)
        return v1 + v2, gv1 + gv2


class Sub(Binary):
    def __call__(self, coordinates):
        v1, gv1 = self.iic1(coordinates)
        v2, gv2 = self.iic2(coordinates)
        return v1 - v2, gv1 - gv2


class Delta(Binary):
    def __call__(self, coordinates):
        b, gb = self.iic1(coordinates)
        e, ge = self.iic2(coordinates)
        return e - b, ge - gb


class Mul(Binary):
    def __call__(self, coordinates):
        v1, gv1 = self.iic1(coordinates)
        v2, gv2 = self.iic2(coordinates)
        return v1*v2, v2*gv1+v1*gv2
        

class Div(Binary):
    def __call__(self, coordinates):
        v1, gv1 = self.iic1(coordinates)
        v2, gv2 = self.iic2(coordinates)
        return v1/v2, (v2*gv1+v1*gv2)/(v2*v2)


class Dot(Binary):
    def __call__(self, coordinates):
        e1, ge1 = self.iic1(coordinates)
        e2, ge2 = self.iic2(coordinates)
        dot = Numeric.dot(e1, e2)
        gdot = Numeric.zeros(coordinates.shape, Numeric.Float)
        for i in range(3):
            gdot[:,i] = e2[i]*ge1[:,i] + e1[i]*ge2[:,i]
        return dot,gdot


class Cos(Binary):
    def __init__(self, distance1, distance2, **keyvals):
        assert isinstance(distance1, Distance), "The first iic must be a distance."
        assert isinstance(distance2, Distance), "The second iic must be a distance."
        Binary.__init__(self, distance1, distance2, **keyvals)
        
    def __call__(self, coordinates):
        d1, gd1 = self.iic1(coordinates)
        d2, gd2 = self.iic2(coordinates)
        e1, ge1 = self.iic1.iic(coordinates)
        e2, ge2 = self.iic2.iic(coordinates)
        edot = Numeric.dot(e1, e2)
        dprod = (d1*d2)
        ds1 = d1*d1
        ds2 = d2*d2
        cos = edot/dprod
        gcos = Numeric.zeros(coordinates.shape, Numeric.Float)
        for i in range(3):
            gcos[:,i] = (e2[i] - edot*e1[i]/ds1)*ge1[:,i] + (e1[i] - edot*e2[i]/ds2)*ge2[:,i]
        gcos /= dprod
        return cos,gcos        

        
class Unary(InternalCoordinate):
    """
    The base class for internal coordinates that take _one_ internal
    coordinate as input parameter. 'iic' in the code stands for "input 
    internal coordinate".
    """

    def __init__(self, iic, **keyvals):
        InternalCoordinate.__init__(self, **keyvals)
        self.iic = iic
        
    def label(self):
        return "%s(%s)" % (self.symbol(), self.iic.label())


class Distance(Unary):
    """
    The Distance ic requires a input internal coordinate that follows
    convetion B. (See class InternalCoordinate for more information about
    these conventions.)
    """

    def __call__(self, coordinates):
        e, ge = self.iic(coordinates)
        distance = math.sqrt(Numeric.dot(e, e))
        gdistance = copy.deepcopy(ge)
        for i in range(3):
            gdistance[:,i] *= e[i]/distance
        return distance, gdistance


class DistanceSqr(Unary):
    """
    The DistanceSqr ic requires a input internal coordinate that follows
    convetion B. (See class InternalCoordinate for more information about
    these conventions.)
    """

    def __call__(self, coordinates):
        e, ge = self.iic(coordinates)
        distancesqr = Numeric.dot(e, e)
        gdistancesqr = copy.deepcopy(ge)
        for i in range(3):
            gdistancesqr[:,i] *= 2*e[i]
        return distancesqr, gdistancesqr


class Sqr(Unary):
    def __call__(self, coordinates):
        e, ge = self.iic(coordinates)
        sqr = e*e
        gsqr = copy.deepcopy(ge)
        gsqr *= 2*e
        return sqr, gsqr


class Sqrt(Unary):
    def __call__(self, coordinates):
        e, ge = self.iic(coordinates)
        sqrt = math.sqrt(e)
        gsqrt = copy.deepcopy(ge)
        gsqrt /= 2*sqrt
        return sqrt, gsqrt


class ArcCos(Unary):
    def __call__(self, coordinates):
        cos, gcos = self.iic(self, coordinates)
        return math.acos(cos), -1 / math.sqrt(1 - x*x)


# Tools for dealing with large sets of internal coordinates


class Collection(object):
    """
    Collection has a twofold goal: (i) Ease the mass creation of internal 
    coordinates and (ii) make sure an internal coordinate is only created once.
    """
    def __init__(self, molecule):
        self.internal_coordinates = {}
        self.molecular_graph = MolecularGraph(molecule)
        self.user_coordinates = {}
        
    def __len__(self):
        return len(self.user_coordinates)
        
    def __getitem__(self, key):
        return self.user_coordinates.get(key)
        
    def __iter__(self):
        return iter(self.user_coordinates)
        
    def add(self, InternalCoordinateClass, *parameters, **keyvals):
        """
        Creates an internal coordinate and adds it to the collection.
        
        This method assures that all the parameters are internal coordinates 
        that are availlable in the collection and the the given internal
        coordinate is not recreated if it already exists.
        """
        for parameter in parameters:
            if isinstance(parameter, InternalCoordinate):
                assert parameter.label() in self.internal_coordinates
        
        internal_coordinate = InternalCoordinateClass(*parameters, **keyvals)
        label = internal_coordinate.label()
        existing_internal_coordinate = self.internal_coordinates.get(label)
        if existing_internal_coordinate == None:
            self.internal_coordinates[internal_coordinate.label()] = internal_coordinate
            return internal_coordinate
        else:
            return existing_internal_coordinate
    
    def add_internal_coordinate(self, tag, internal_coordinate):
        internal_coordinate.tag = tag
        existing_internal_coordinates = self.user_coordinates.get(tag)
        if existing_internal_coordinates == None:
            self.user_coordinates[tag] = [internal_coordinate]
        else:
            existing_internal_coordinates.append(internal_coordinate)
   
    def add_long_range_distances(self, criteria_sets):
        def all_pairs(atom_criteria):
            first_criterium = atom_criteria[0]
            first_criterium.set_molecular_graph(self.molecular_graph)
            second_criterium = atom_criteria[1]
            second_criterium.set_molecular_graph(self.molecular_graph)
            molecule = self.molecular_graph.molecule
            for index1 in xrange(len(molecule.numbers)):
                for index2 in xrange(index1+1, len(molecule.numbers)):
                    if first_criterium(index1) and second_criterium(index2):
                        yield (index1, index2)
                    elif first_criterium(index2) and second_criterium(index1):
                        yield (index2, index1)
                
    
        nonbonded_pairs = dict((tag, set(all_pairs(atom_criteria))) for tag, atom_criteria, bond_criteria, filter_tags in criteria_sets.yield_criteria())

        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets.bond_excludes()):
            id = tuple([match.get_destination(source) for source in [0, 1]])
            nonbonded_pairs[tag].discard(id)
            reverse_id = (id[1], id[0])
            nonbonded_pairs[tag].discard(reverse_id)

        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets.bend_excludes()):
            id = tuple([match.get_destination(source) for source in [0, 2]])
            nonbonded_pairs[tag].discard(id)
            reverse_id = (id[1], id[0])
            nonbonded_pairs[tag].discard(reverse_id)

        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets.dihedral_excludes()):
            id = tuple([match.get_destination(source) for source in [0, 3]])
            nonbonded_pairs[tag].discard(id)
            reverse_id = (id[1], id[0])
            nonbonded_pairs[tag].discard(reverse_id)

        def distance(id, tag):
            s0 = self.add(Select, id[0])
            s1 = self.add(Select, id[1])
            e = self.add(Delta, s1, s0)
            d = self.add(
                Distance, 
                e, 
                name="long range distance %i-%i (%s)" % (id + (tag,)),
                id=id
            )
            return d

        for tag, ids in nonbonded_pairs.iteritems():
            for id in ids:
                self.add_internal_coordinate(tag, distance(id, tag))
   
    def add_bond_lengths(self, criteria_sets):
        """
        Adds the bond lengths described in criteria_sets to the collection.
        
        Arguments
        criteria_sets -- see pychem.molecular_graphs
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1]])
            s0 = self.add(Select, id[0])
            s1 = self.add(Select, id[1])
            e = self.add(Delta, s1, s0)
            d = self.add(
                Distance, 
                e, 
                name="bond length %i-%i (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, d)
            
    def add_bend_cosines(self, criteria_sets):
        """
        Adds the cosines of the bend angles described in criteria_sets to the 
        collection.
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1, 2]])
            s0 = self.add(Select, id[0])
            s1 = self.add(Select, id[1])
            s2 = self.add(Select, id[2])
            e1 = self.add(Delta, s1, s0)
            e2 = self.add(Delta, s1, s2)
            d1 = self.add(Distance, e1)
            d2 = self.add(Distance, e2)
            c = self.add(
                Cos, 
                d1,
                d2, 
                name="bend cos %i-%i-%i (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, c)
        
    def add_bend_spans(self, criteria_sets):
        """
        Adds the distances that span the bend angles described in criteria_sets
        to the collection.
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1, 2]])
            s0 = self.add(Select, id[0])
            s2 = self.add(Select, id[2])
            e = self.add(Delta, s0, s2)
            d = self.add(
                Distance, 
                e, 
                name="bend span %i-%i-%i (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, d)

    def add_dihedral_cosines(self, criteria_sets):
        """
        Adds the dihedral angles described in criteria_sets to the collection.
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1, 2, 3]])
            s0 = self.add(Select, id[0])
            s1 = self.add(Select, id[1])
            s2 = self.add(Select, id[2])
            s3 = self.add(Select, id[3])
            el = self.add(Delta, s1, s0)
            em = self.add(Delta, s1, s2)
            er = self.add(Delta, s2, s3)
            dll = self.add(DistanceSqr, el)
            dmm = self.add(DistanceSqr, em)
            drr = self.add(DistanceSqr, er)
            dlm = self.add(Dot, el, em)
            dmr = self.add(Dot, em, er)
            drl = self.add(Dot, er, el)
            #
            tin = self.add(Mul, dlm, dmr)
            tout = self.add(Mul, drl, dmm)
            t = self.add(Sub, tout, tin)
            noutl = self.add(Mul, dll, dmm)
            noutr = self.add(Mul, dmm, drr)
            dlm2 = self.add(Sqr, dlm)
            dmr2 = self.add(Sqr, dmr)
            nl = self.add(Sub, noutl, dlm2)
            nr = self.add(Sub, noutr, dmr2)
            n2 = self.add(Mul, nl, nr)
            n = self.add(Sqrt, n2)
            dihedral_cos = self.add(
                Div, 
                t, 
                n,
                name="dihedral cos %i-%i-%i-%i (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, dihedral_cos)

    def add_dihedral_spans(self, criteria_sets):
        """
        Adds the distances that span the dihedral angles described in
        criteria_sets to the collection.
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1, 2, 3]])
            s0 = self.add(Select, id[0])
            s3 = self.add(Select, id[3])
            e = self.add(Delta, s0, s3)
            d = self.add(
                Distance, 
                e, 
                name="dihedral span %i-%i-%i-%i (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, d)

    def add_out_of_plane_cosines(self, criteria_sets):
        """
        Adds the out of plane cosines described in criteria_sets to the
        collection.
        """
        result = dict((tag, []) for tag in criteria_sets.yield_tags())
        for tag, match in self.molecular_graph.yield_subgraphs(criteria_sets):
            id = tuple([match.get_destination(source) for source in [0, 1, 2, 3]])
            s0 = self.add(Select, id[0])
            s1 = self.add(Select, id[1])
            s2 = self.add(Select, id[2])
            s3 = self.add(Select, id[3])
            rt = self.add(Delta, s1, s0)
            ra = self.add(Delta, s1, s2)
            rb = self.add(Delta, s1, s3)
            dat = self.add(Dot, ra, rt)
            dbt = self.add(Dot, rb, rt)
            dl = self.add(Mul, rb, dat)
            dr = self.add(Mul, ra, dbt)
            d = self.add(Sub, dl, dr)
            ddb = self.add(Dot, d, rb)
            dda = self.add(Dot, d, ra)
            pl = self.add(Mul, ddb, ra)
            pr = self.add(Mul, dda, rb)
            p = self.add(Sub, pl, pr)
            t = self.add(Dot, p, rt)
            n1 = self.add(DistanceSqr, p)
            n2 = self.add(DistanceSqr, rt)
            nn = self.add(Mul, n1, n2)
            n = self.add(Sqrt, nn)
            out_of_plane_cos = self.add(
                Div,
                t,
                n,
                name="out of plane cos %i-%i (%i,%i) (%s)" % (id + (tag,)),
                id=id
            )
            self.add_internal_coordinate(tag, out_of_plane_cos) 

    def yield_related_internal_coordinates(self, tag1, tag2, order_related=2):
        for ic1 in self[tag1]:
            for ic2 in self[tag2]:
                if len(set(ic1.id) & set(ic2.id)) >= order_related:
                    yield (ic1, ic2)

    def add_related_products(self, tag, tag1, tag2, order_related=2):
        result = []
        for ic1, ic2 in self.yield_related_internal_coordinates(tag1, tag2, order_related):
            product = self.add(
                Mul,
                ic1,
                ic2,
                name="%s x %s (%s)" % (ic1.name, ic2.name, tag)
            )
            self.add_internal_coordinate(tag, product)


# Tools for solving the jacobian system

class Configuration(object):
    pass


class JacobianSolver(object):
    """
    The JacobianSolver is a utility class solves the jacobian system and saves
    additional calculations about a single point calculation into a
    Configuration object.
    """
    
    def __init__(self, internal_coordinates, num_atoms):
        assert len(internal_coordinates) >= num_atoms*3-6, "Specify more internal coordinates"
        self.internal_coordinates = internal_coordinates
        # internal coordinates mask: This masks out the non-independant
        # carthesian coordinates in the center of mass frame. It is assumed
        # that the molecule is transformed into it's normalized frame.
        # (see ffaudit.molecules.Molecule.normalize())
        # the masked coordinates are: x1, y1, z1, y2, z2, z3
        self.internal_mask = Numeric.ones((num_atoms, 3))
        self.internal_mask[0,:] = 0
        self.internal_mask[1,1:] = 0
        self.internal_mask[2,2:] = 0
        self.internal_mask = Numeric.ravel(self.internal_mask)
        # One must specify at least as much internal coordinates as
        # there are internal degrees of freedom in a molecule.
        
    def analyze_gradient(self, job):
        """
        Extract all information about partial derivates towards new internal 
        coordinates of a single point job.
        
        The fields in the configuration object are interpreted as:
        ai_energy -- the ab initio energy of the sp job
        values -- the values of the internal coordinates
        jacobian -- the masked jacobian (without redundant carthesian
                    coordinates)
        full_jacobian -- with redundant carthesian coordinatess, and with a
                         different shape: rang-2 (num_ic x 3*num_carth)
        V, S, Wt -- see Singular value decomposition (for solving general
                    linear systems)
        W -- the transpose of Wt
        rank -- the number of independent internal coordinates
        gradient -- the masked carthesian gradient
        full_gradient -- the unmasked carthesian gradient rang-1 (num_cart*3)
        particular -- a particular solution for the partial derivates towards
                      a set of redundant internal coordinates (a solution of the
                      jacobian system)
        null_space -- the null_space of the jacobian system
        """
        result = Configuration()
        result.ai_energy = job.energy

        jacobian = []
        full_jacobian = []
        values = []
        for internal_coordinate in self.internal_coordinates:
            value, pd = internal_coordinate(job.input_molecule.coordinates)
            values.append(value)
            jacobian.append(Numeric.compress(self.internal_mask, Numeric.ravel(pd)))
            full_jacobian.append(Numeric.ravel(pd))
            
        result.values = Numeric.array(values)
        result.jacobian = Numeric.transpose(Numeric.array(jacobian))
        result.full_jacobian = Numeric.transpose(Numeric.array(full_jacobian))
        
        result.V, result.S, result.Wt = LinearAlgebra.singular_value_decomposition(result.jacobian, True)
        result.W = Numeric.transpose(result.Wt)
        result.rank = len(result.S)
        
        result.gradient = Numeric.compress(self.internal_mask, Numeric.ravel(job.gradient))
        result.full_gradient = Numeric.ravel(job.gradient)
        result.particular = Numeric.dot(result.W[:,:result.rank], Numeric.transpose(Numeric.transpose(Numeric.dot(result.gradient, result.V))/result.S))

        if result.W.shape[1] > result.rank:
            result.nullspace = result.W[:,result.rank:]
            # initial guess = solution with minimum norm
            result.particular -= Numeric.dot(result.nullspace, Numeric.dot(result.particular, result.nullspace))
            
        return result


    def analyze_hessian(self, job, molecule):
        result = Configuration()
        full_jacobian = []
        for internal_coordinate in self.internal_coordinates:
            value, pd = internal_coordinate(molecule.coordinates)
            full_jacobian.append(Numeric.ravel(pd))
            
        full_jacobian = Numeric.transpose(Numeric.array(full_jacobian))
        normalized_full_jacobian = full_jacobian.copy()
        for col in Numeric.transpose(normalized_full_jacobian):
            col[:] /= math.sqrt(Numeric.dot(col, Numeric.dot(job.hessian, col)))
        
        coupling = Numeric.dot(Numeric.transpose(normalized_full_jacobian), Numeric.dot(job.hessian, normalized_full_jacobian))
        return coupling
        
