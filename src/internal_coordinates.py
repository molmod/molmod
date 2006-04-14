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
import numpy

__all__ = [
    "InternalCoordinate", "Select", 
    "Binary", "Add", "Sub", "Delta", "Mul", "Div", "Dot",
    "Unary", "Distance", "DistanceSqr", "Sqr", "Sqrt",
    "InternalCoordinatesCache"
]


# Classes for calculating the partial derivatives of the internal
# coordinates towards the cartesian coordinates

SCALAR = 0
VECTOR = 2


class InternalCoordinate(object):
    def __init__(self, output_style, **keyvals):
        self.__dict__.update(keyvals)
        self.output_style = output_style
        
    def description(self):
        temp = str(self.__class__)
        return temp[temp.rfind(".")+1:-2].lower()
        
    def label(self):
        raise NotImplementedError

    def value_tangent(self, coordinates):
        raise NotImplementedError
        # interpretation of tangent:
        #    - In case the returned value is a vector (as a function of the coordinates)
        #         array[i, k, m] = partial derivative of result[m] towards coordinate k of atom i
        #    - In case the returned value is a scalar (...)
        #         array[i, k] = partial derivative of result towards coordinate k of atom i

    def curvature(self, coordinates):
        raise NotImplementedError
        # interpretation of tangent:
        #    - In case the returned value is a vector (as a function of the coordinates)
        #         array[i, k, j, l, m] = second order partial derivative of result[m] towards coordinate k of atom i and coordinate l of atom j
        #    - In case the returned value is a scalar (...)
        #         array[i, k, j, l] = second order partial derivative of result towards coordinate k of atom i and coordinate l of atom j


class Select(InternalCoordinate):
    """
    This internal coordinate is only a 'helper' for calculating other
    internal coordinates. It selects a given atom from the molecule.
    """

    def __init__(self, index, **keyvals):
        InternalCoordinate.__init__(self, VECTOR, **keyvals)
        self.index = index
        
    def value_tangent(self, coordinates):
        tangent = numpy.zeros((len(coordinates), 3, 3), float)
        tangent[self.index, 0, 0] = 1
        tangent[self.index, 1, 1] = 1
        tangent[self.index, 2, 2] = 1
        return coordinates[self.index], tangent

    def curvature(self, coordinates):
        curvature = numpy.zeros((len(coordinates), 3, len(coordinates), 3, 3), float)
        return curvature
        
    def label(self):
        return "%s(%i)" % (self.description(), self.index)


class Binary(InternalCoordinate):
    """
    The base class for internal coordinates that take _two_ internal
    coordinate as input parameter. 'iic' in the code stands for "input 
    internal coordinate".
    """
       
    def __init__(self, output_style, iic1, iic2, **keyvals):
        InternalCoordinate.__init__(self, output_style, **keyvals)
        self.iic1 = iic1
        self.iic2 = iic2
        
    def label(self):
        return "%s(%s,%s)" % (self.description(), self.iic1.label(), self.iic2.label())


class HybridBinary(Binary):
    def __init__(self, iic1, iic2, **keyvals):
        assert iic1.output_style == iic2.output_style
        Binary.__init__(self, iic1.output_style, iic1, iic2, **keyvals)


class Add(HybridBinary):
    def value_tangent(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        return v1 + v2, t1 + t2

    def curvature(self, coordinates):
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        return c1 + c2


class Sub(HybridBinary):
    def value_tangent(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        return v1 - v2, t1 - t2

    def curvature(self, coordinates):
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        return c1 - c2


class Delta(HybridBinary):
    def value_tangent(self, coordinates):
        b, tb = self.iic1.value_tangent(coordinates)
        e, te = self.iic2.value_tangent(coordinates)
        return e - b, te - tb

    def curvature(self, coordinates):
        b = self.iic1.curvature(coordinates)
        e = self.iic2.curvature(coordinates)
        return e - b


class ScalarBinary(Binary):
    def __init__(self, iic1, iic2, **keyvals):
        assert iic1.output_style == SCALAR
        assert iic2.output_style == SCALAR
        Binary.__init__(self, SCALAR, iic1, iic2, **keyvals)


class Mul(ScalarBinary):
    def value_tangent(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        return v1*v2, v2*t1+v1*t2

    def curvature(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        return v1*c2 + v2*c1 + numpy.multiply.outer(t1, t2) + numpy.multiply.outer(t2, t1)
        

class Div(ScalarBinary):
    def value_tangent(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        return v1/v2, (v2*t1-v1*t2)/(v2*v2)

    def curvature(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        v2, t2 = self.iic2.value_tangent(coordinates)
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        return (c1 - (
            numpy.multiply.outer(t1, t2)
            + numpy.multiply.outer(t2, t1)
            + v1*c2
            - 2*v1*numpy.multiply.outer(t2, t2)/v2
        )/v2)/v2
        

class Dot(Binary):
    def __init__(self, iic1, iic2, **keyvals):
        assert iic1.output_style == VECTOR
        assert iic2.output_style == VECTOR
        Binary.__init__(self, SCALAR, iic1, iic2, **keyvals)

    def value_tangent(self, coordinates):
        e1, t1 = self.iic1.value_tangent(coordinates)
        e2, t2 = self.iic2.value_tangent(coordinates)
        dot = numpy.dot(e1, e2)
        tdot = numpy.zeros(coordinates.shape, float)
        for i in range(3):
            tdot += e2[i]*t1[:,:,i] + e1[i]*t2[:,:,i]
        return dot, tdot

    def curvature(self, coordinates):
        e1, t1 = self.iic1.value_tangent(coordinates)
        e2, t2 = self.iic2.value_tangent(coordinates)
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        cdot = numpy.zeros(coordinates.shape + coordinates.shape, float)
        for i in range(3):
            cdot += e1[i]*c2[:,:,:,:,i] + numpy.multiply.outer(t1[:,:,i], t2[:,:,i]) + numpy.multiply.outer(t2[:,:,i], t1[:,:,i]) + e2[i]*c1[:,:,:,:,i]
        return cdot
        

class Scale(Binary):
    def __init__(self, iic1, iic2, **keyvals):
        assert iic1.output_style == SCALAR
        assert iic2.output_style == VECTOR
        Binary.__init__(self, VECTOR, iic1, iic2, **keyvals)

    def value_tangent(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        e2, t2 = self.iic2.value_tangent(coordinates)
        tscale = v1*t2
        for i in xrange(3):
            tscale[:,:,i] += e2[i]*t1
        return v1*e2, tscale
        
    def curvature(self, coordinates):
        v1, t1 = self.iic1.value_tangent(coordinates)
        e2, t2 = self.iic2.value_tangent(coordinates)
        c1 = self.iic1.curvature(coordinates)
        c2 = self.iic2.curvature(coordinates)
        cscale = v1*c2
        for i in range(3):
            cscale[:,:,:,:,i] = c1*e2[i] + numpy.multiply.outer(t1, t2[:,:,i]) + numpy.multiply.outer(t2[:,:,i], t1)
        return cscale


class Unary(InternalCoordinate):
    """
    The base class for internal coordinates that take _one_ internal
    coordinate as input parameter. 'iic' in the code stands for "input 
    internal coordinate".
    """

    def __init__(self, output_style, iic, **keyvals):
        InternalCoordinate.__init__(self, output_style, **keyvals)
        self.iic = iic
        
    def label(self):
        return "%s(%s)" % (self.description(), self.iic.label())


class Factor(Unary):
    def __init__(self, factor, iic, **keyvals):
        Unary.__init__(self, iic.output_style, iic, **keyvals)
        self.factor = factor
        
    def value_tangent(self, coordinates):
        e, t = self.iic.value_tangent(coordinates)
        return self.factor*e, self.factor*t
        
    def curvature(self, coordinates):
        return self.factor*self.iic.curvature(coordinates)

    def label(self):
        return "(%s*%s)" % (self.factor, self.iic.label())


class Measure(Unary):
    def __init__(self, iic, **keyvals):
        assert iic.output_style == VECTOR
        Unary.__init__(self, SCALAR, iic, **keyvals)


class Distance(Measure):
    """
    The Distance ic requires a input internal coordinate that follows
    convetion B. (See class InternalCoordinate for more information about
    these conventions.)
    """

    def value_tangent(self, coordinates):
        e, te = self.iic.value_tangent(coordinates)
        distance = math.sqrt(numpy.dot(e, e))
        tdistance = numpy.zeros(coordinates.shape, float)
        for i in range(3):
            tdistance += te[:,:,i]*e[i]/distance
        return distance, tdistance

    def curvature(self, coordinates):
        e, t = self.iic.value_tangent(coordinates)
        c = self.iic.curvature(coordinates)
        d2 = numpy.dot(e, e)
        d = math.sqrt(d2)
        cd = numpy.zeros(coordinates.shape + coordinates.shape, float)
        for i in range(3):
            cd += (numpy.multiply.outer(t[:,:,i], t[:,:,i]) + e[i]*c[:,:,:,:,i])/d
            for j in range(3):
                cd -= (numpy.multiply.outer(t[:,:,i], t[:,:,j])*e[i]*e[j])/(d*d2)
        return cd


class DistanceSqr(Measure):
    """
    The DistanceSqr ic requires a input internal coordinate that follows
    convetion B. (See class InternalCoordinate for more information about
    these conventions.)
    """

    def value_tangent(self, coordinates):
        e, te = self.iic.value_tangent(coordinates)
        distancesqr = numpy.dot(e, e)
        tdistancesqr = numpy.zeros(coordinates.shape, float)
        for i in range(3):
            tdistancesqr += 2*te[:,:,i]*e[i]
        return distancesqr, tdistancesqr

    def curvature(self, coordinates):
        e, t = self.iic.value_tangent(coordinates)
        c = self.iic.curvature(coordinates)
        cd2 = numpy.zeros(coordinates.shape + coordinates.shape, float)
        for i in range(3):
            cd2 += 2*(numpy.multiply.outer(t[:,:,i], t[:,:,i]) + e[i]*c[:,:,:,:,i])
        return cd2


class ScalarUnary(Unary):
    def __init__(self, iic, **keyvals):
        assert iic.output_style == SCALAR
        Unary.__init__(self, SCALAR, iic, **keyvals)


class Sqr(ScalarUnary):
    def value_tangent(self, coordinates):
        v, t = self.iic.value_tangent(coordinates)
        return v*v, 2*v*t

    def curvature(self, coordinates):
        v, t = self.iic.value_tangent(coordinates)
        c = self.iic.curvature(coordinates)
        return 2*(numpy.multiply.outer(t, t) + v*c)


class Sqrt(ScalarUnary):
    def value_tangent(self, coordinates):
        v, t = self.iic.value_tangent(coordinates)
        sqrt = math.sqrt(v)
        return sqrt, t/(2*sqrt)

    def curvature(self, coordinates):
        v, t = self.iic.value_tangent(coordinates)
        sqrt = math.sqrt(v)
        c = self.iic.curvature(coordinates)
        return 0.5*c/sqrt - 0.25*numpy.multiply.outer(t, t)/(v*sqrt)


# Tools for dealing with large sets of internal coordinates


class InternalCoordinatesCache(object):
    """
    InternalCoordinatesCache has a twofold goal: (i) Ease the mass creation of 
    internal coordinates and (ii) make sure an internal coordinate is only
    created once.
    """
    def __init__(self, molecule):
        self.internal_coordinates = {}
        self.molecular_graph = MolecularGraph(molecule)
        self.user_coordinates = {}
        
    def __getitem__(self, key):
        return self.user_coordinates[key]
        
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
                name="long range distance %i-%i" % id,
                id=id,
                symbol="D%i.%i" % id
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
                name="bond length %i-%i" % id,
                id=id,
                symbol="D%i-%i" % id
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
            dot = self.add(Dot, e1, e2)
            dd = self.add(Mul, d1, d2)
            c = self.add(
                Div, 
                dot,
                dd, 
                name="bend cos %i-%i-%i" % id,
                id=id,
                symbol="C%i-%i-%i" % id
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
                name="bend span %i-%i-%i" % id,
                id=id,
                symbol="D%i^%i" % (id[0], id[2])
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
                name="dihedral cos %i-%i-%i-%i" % id,
                id=id,
                symbol="C%i-%i-%i-%i" % id
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
                name="dihedral span %i-%i-%i-%i" % id,
                id=id,
                symbol="D%i~%i" % (id[0], id[3])
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
            dl = self.add(Scale, dat, rb)
            dr = self.add(Scale, dbt, ra)
            d = self.add(Sub, dl, dr)
            ddb = self.add(Dot, d, rb)
            dda = self.add(Dot, d, ra)
            pl = self.add(Scale, ddb, ra)
            pr = self.add(Scale, dda, rb)
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
                name="out of plane cos %i-%i (%i,%i)" % id,
                id=id,
                symbol="C%i-%i(%i,%i)" % id
            )
            self.add_internal_coordinate(tag, out_of_plane_cos) 


    def add_out_of_plane_distances(self, criteria_sets):
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
            ra = self.add(Delta, s1, s2)
            rb = self.add(Delta, s1, s3)
            rc = self.add(Delta, s1, s0)
            daa = self.add(DistanceSqr, ra)
            dbb = self.add(DistanceSqr, rb)
            dcc = self.add(DistanceSqr, rc)
            dab = self.add(Dot, ra, rb)
            dac = self.add(Dot, ra, rc)
            dbc = self.add(Dot, rb, rc)
            paabb = self.add(Mul, daa, dbb)
            pab2 = self.add(Sqr, dab)
            dnn = self.add(Sub, paabb, pab2)
            rt1 = self.add(Scale, dac, rb)
            rt2 = self.add(Scale, dbc, ra)
            re = self.add(Sub, rt1, rt2)
            dee = self.add(DistanceSqr, re)
            t3 = self.add(Div, dee, dnn)
            dout2 = self.add(Sub, dcc, t3)
            out_of_plane_distance = self.add(
                Sqrt,
                dout2,
                name="out of plane distance %i-(%i,%i,%i)" % id,
                id=id,
                symbol="D%i-(%i,%i,%i)" % id
            )
            self.add_internal_coordinate(tag, out_of_plane_distance) 

    def yield_related_internal_coordinates(self, tag1, tag2, order_related=2):
        for ic1 in self[tag1]:
            for ic2 in self[tag2]:
                if tag1==tag2 and ic1 >= ic2:
                    continue
                num_common = len(set(ic1.id) & set(ic2.id))
                if num_common >= order_related:
                    yield (ic1, ic2)

    def add_related_products(self, tag, tag1, tag2, order_related=2):
        result = []
        for ic1, ic2 in self.yield_related_internal_coordinates(tag1, tag2, order_related):
            product = self.add(
                Mul,
                ic1,
                ic2,
                name="%s x %s" % (ic1.name, ic2.name),
                symbol="%s*%s" % (ic1.symbol, ic2.symbol)
            )
            self.add_internal_coordinate(tag, product)
